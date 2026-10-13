// SPDX-License-Identifier: GPL-2.0-or-later

use crate::grid::Grid2D;

/// Forme d'onde d'une source ponctuelle.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Waveform {
    /// Ondelette de Ricker (dérivée seconde de gaussienne), centrée à `t = 1/freq`.
    /// Sans composante continue : idéale pour la simulation.
    Ricker { freq: f32 },
    /// Sinusoïde continue.
    Sine { freq: f32 },
}

impl Waveform {
    /// Valeur du signal au temps `t` (secondes).
    pub fn eval(&self, t: f64) -> f32 {
        match *self {
            Waveform::Ricker { freq } => {
                let f = freq as f64;
                let a = (std::f64::consts::PI * f * (t - 1.0 / f)).powi(2);
                ((1.0 - 2.0 * a) * (-a).exp()) as f32
            }
            Waveform::Sine { freq } => (2.0 * std::f64::consts::PI * freq as f64 * t).sin() as f32,
        }
    }
}

#[derive(Debug, Clone)]
struct Source {
    x: usize,
    y: usize,
    waveform: Waveform,
    amplitude: f32,
}

/// Solveur acoustique 2D (FDTD, schéma saute-mouton sur grille décalée).
///
/// Disposition : `p[x,y]` au centre des cellules, `vx[x,y]` entre `p[x,y]` et
/// `p[x+1,y]`, `vy[x,y]` entre `p[x,y]` et `p[x,y+1]`.
/// Les frontières sont des murs rigides (vitesse normale nulle), sauf si
/// [`set_absorbing_boundary`](Self::set_absorbing_boundary) est utilisé.
pub struct AcousticSolver2D {
    pub p: Grid2D,
    pub vx: Grid2D,
    pub vy: Grid2D,
    rho: Grid2D,
    inv_rho: Grid2D,
    bulk_modulus: Grid2D,
    /// Facteur multiplicatif appliqué à chaque pas (1.0 = pas d'amortissement).
    damp: Grid2D,
    // Nouveaux champs pour les PML
    sigma_x: Grid2D,
    sigma_y: Grid2D,

    zero_row: Vec<f32>,
    sources: Vec<Source>,
    receivers: Vec<(usize, usize)>,
    traces: Vec<Vec<f32>>,
    time: f64,
}

impl AcousticSolver2D {
    pub fn new(nx: usize, ny: usize, dx: f32, dy: f32, default_rho: f32, default_c: f32) -> Self {
        assert!(nx >= 2 && ny >= 2, "la grille doit faire au moins 2x2");
        assert!(
            default_rho > 0.0 && default_c > 0.0,
            "rho et c doivent être > 0"
        );

        let mut rho = Grid2D::new(nx, ny, dx, dy);
        rho.fill(default_rho);
        let mut inv_rho = Grid2D::new(nx, ny, dx, dy);
        inv_rho.fill(1.0 / default_rho);
        let mut bulk_modulus = Grid2D::new(nx, ny, dx, dy);
        bulk_modulus.fill(default_rho * default_c * default_c);
        let mut damp = Grid2D::new(nx, ny, dx, dy);
        damp.fill(1.0);

        // Initialisation des grilles PML à zéro par défaut
        let sigma_x = Grid2D::new(nx, ny, dx, dy);
        let sigma_y = Grid2D::new(nx, ny, dx, dy);

        Self {
            p: Grid2D::new(nx, ny, dx, dy),
            vx: Grid2D::new(nx, ny, dx, dy),
            vy: Grid2D::new(nx, ny, dx, dy),
            rho,
            inv_rho,
            bulk_modulus,
            damp,
            sigma_x, // Ajouté ici
            sigma_y, // Ajouté ici
            zero_row: vec![0.0; nx],
            sources: Vec::new(),
            receivers: Vec::new(),
            traces: Vec::new(),
            time: 0.0,
        }
    }

    // ---------------------------------------------------------------- milieu

    /// Densité (lecture seule). Utiliser [`set_medium`](Self::set_medium) pour modifier.
    pub fn rho(&self) -> &Grid2D {
        &self.rho
    }

    /// Configure des couches absorbantes PML sur les bords.
    /// `width` : épaisseur en nombre de cellules (ex: 20 à 30).
    /// `max_sigma` : intensité maximale de l'amortissement.
    pub fn set_pml_boundary(&mut self, width: usize, max_sigma: f32) {
        let (nx, ny) = (self.p.nx, self.p.ny);

        for y in 0..ny {
            for x in 0..nx {
                // Distance par rapport aux bords pour l'axe X
                let dist_x_left = if x < width { width - x } else { 0 };
                let dist_x_right = if x >= nx - width {
                    x - (nx - 1 - width)
                } else {
                    0
                };
                let dist_x = dist_x_left.max(dist_x_right) as f32;

                if width > 0 && dist_x > 0.0 {
                    let t = dist_x / width as f32;
                    let val = max_sigma * t * t; // Profil quadratique
                    self.sigma_x.set(x, y, val);
                } else {
                    self.sigma_x.set(x, y, 0.0);
                }

                // Distance par rapport aux bords pour l'axe Y
                let dist_y_top = if y < width { width - y } else { 0 };
                let dist_y_bottom = if y >= ny - width {
                    y - (ny - 1 - width)
                } else {
                    0
                };
                let dist_y = dist_y_top.max(dist_y_bottom) as f32;

                if width > 0 && dist_y > 0.0 {
                    let t = dist_y / width as f32;
                    let val = max_sigma * t * t;
                    self.sigma_y.set(x, y, val);
                } else {
                    self.sigma_y.set(x, y, 0.0);
                }
            }
        }
    }

    /// Module de compressibilité K = ρc² (lecture seule).
    pub fn bulk_modulus(&self) -> &Grid2D {
        &self.bulk_modulus
    }

    /// Définit le milieu d'une cellule : densité `rho` et vitesse du son `c`.
    pub fn set_medium(&mut self, x: usize, y: usize, rho: f32, c: f32) {
        assert!(rho > 0.0 && c > 0.0, "rho et c doivent être > 0");
        self.rho.set(x, y, rho);
        self.inv_rho.set(x, y, 1.0 / rho);
        self.bulk_modulus.set(x, y, rho * c * c);
    }

    /// Définit le milieu d'un rectangle `[x0, x1) × [y0, y1)`.
    pub fn set_medium_rect(
        &mut self,
        x0: usize,
        y0: usize,
        x1: usize,
        y1: usize,
        rho: f32,
        c: f32,
    ) {
        for y in y0..y1.min(self.p.ny) {
            for x in x0..x1.min(self.p.nx) {
                self.set_medium(x, y, rho, c);
            }
        }
    }

    /// Vitesse du son maximale dans le domaine.
    pub fn c_max(&self) -> f32 {
        self.bulk_modulus
            .as_slice()
            .iter()
            .zip(self.inv_rho.as_slice())
            .map(|(k, ir)| (k * ir).sqrt())
            .fold(0.0, f32::max)
    }

    /// Pas de temps maximal stable (condition CFL 2D) pour le milieu actuel.
    pub fn max_stable_dt(&self) -> f32 {
        let inv = (1.0 / (self.p.dx * self.p.dx) + 1.0 / (self.p.dy * self.p.dy)).sqrt();
        1.0 / (self.c_max() * inv)
    }

    /// `safety × max_stable_dt()`, avec `0 < safety <= 1` (0.9 est un bon choix).
    pub fn stable_dt(&self, safety: f32) -> f32 {
        safety * self.max_stable_dt()
    }

    // ------------------------------------------------------ bords absorbants

    /// Couche éponge de `width` cellules sur les 4 bords.
    ///
    /// `max_damping` (entre 0 et 1) est la fraction d'amplitude retirée **à chaque
    /// pas** sur le bord extérieur ; le profil décroît de façon quadratique vers
    /// l'intérieur. Valeurs typiques : `width` 15–30, `max_damping` 0.05–0.2.
    /// `width = 0` désactive l'amortissement (murs rigides).
    pub fn set_absorbing_boundary(&mut self, width: usize, max_damping: f32) {
        assert!(
            (0.0..1.0).contains(&max_damping),
            "max_damping doit être dans [0, 1)"
        );
        let (nx, ny) = (self.p.nx, self.p.ny);
        for y in 0..ny {
            for x in 0..nx {
                let d = x.min(nx - 1 - x).min(y).min(ny - 1 - y);
                let f = if width > 0 && d < width {
                    let t = (width - d) as f32 / width as f32;
                    1.0 - max_damping * t * t
                } else {
                    1.0
                };
                self.damp.set(x, y, f);
            }
        }
    }

    // -------------------------------------------------- sources, récepteurs

    /// Ajoute une source ponctuelle pilotée par une forme d'onde.
    pub fn add_source(&mut self, x: usize, y: usize, waveform: Waveform, amplitude: f32) {
        assert!(x < self.p.nx && y < self.p.ny, "source hors grille");
        self.sources.push(Source {
            x,
            y,
            waveform,
            amplitude,
        });
    }

    /// Injecte instantanément une impulsion de pression.
    pub fn inject_source(&mut self, x: usize, y: usize, value: f32) {
        let current = self.p.get(x, y);
        self.p.set(x, y, current + value);
    }

    /// Ajoute un récepteur ; retourne son identifiant.
    pub fn add_receiver(&mut self, x: usize, y: usize) -> usize {
        assert!(x < self.p.nx && y < self.p.ny, "récepteur hors grille");
        self.receivers.push((x, y));
        self.traces.push(Vec::new());
        self.receivers.len() - 1
    }

    /// Signal de pression enregistré (un échantillon par pas de temps).
    pub fn receiver(&self, id: usize) -> &[f32] {
        &self.traces[id]
    }

    /// Temps simulé écoulé (secondes).
    pub fn time(&self) -> f64 {
        self.time
    }

    // ------------------------------------------------------------ simulation

    /// Effectue un pas de simulation temporel `dt`.
    pub fn step(&mut self, dt: f32) {
        let nx = self.p.nx;
        let ny = self.p.ny;
        let inv_dx = 1.0 / self.p.dx;
        let inv_dy = 1.0 / self.p.dy;

        // 0. Sources
        for s in &self.sources {
            let idx = self.p.index(s.x, s.y);
            self.p.as_mut_slice()[idx] += s.amplitude * s.waveform.eval(self.time);
        }

        // 1. Vitesses à partir du gradient de pression (avec PML)
        {
            let p = self.p.as_slice();
            let inv_rho = self.inv_rho.as_slice();
            let sig_x = self.sigma_x.as_slice();
            let sig_y = self.sigma_y.as_slice();

            let vel = |y: usize, vxr: &mut [f32], vyr: &mut [f32]| {
                let o = y * nx;
                let prow = &p[o..o + nx];
                let irow = &inv_rho[o..o + nx];
                let sx_row = &sig_x[o..o + nx];
                let vxr = &mut vxr[..nx];

                for i in 0..nx - 1 {
                    let damping_x = 1.0 / (1.0 + 0.5 * dt * sx_row[i]);
                    let dp_dx = (prow[i + 1] - prow[i]) * inv_dx;
                    vxr[i] = (vxr[i] - dt * irow[i] * dp_dx) * damping_x;
                }

                if y + 1 < ny {
                    let pn = &p[o + nx..o + 2 * nx];
                    let sy_next = &sig_y[o + nx..o + 2 * nx];
                    let vyr = &mut vyr[..nx];

                    for i in 0..nx {
                        let damping_y = 1.0 / (1.0 + 0.5 * dt * sy_next[i]);
                        let dp_dy = (pn[i] - prow[i]) * inv_dy;
                        vyr[i] = (vyr[i] - dt * irow[i] * dp_dy) * damping_y;
                    }
                }
            };
            for_each_row2(self.vx.as_mut_slice(), self.vy.as_mut_slice(), nx, vel);
        }

        // 2. Pression à partir de la divergence des vitesses (avec PML)
        {
            let vx = self.vx.as_slice();
            let vy = self.vy.as_slice();
            let bulk = self.bulk_modulus.as_slice();
            let sig_x = self.sigma_x.as_slice();
            let sig_y = self.sigma_y.as_slice();
            let zeros = &self.zero_row[..];

            let prs = |y: usize, prow: &mut [f32]| {
                let o = y * nx;
                let vxr = &vx[o..o + nx];
                let vyr = &vy[o..o + nx];
                let vyp: &[f32] = if y > 0 { &vy[o - nx..o] } else { zeros };
                let kr = &bulk[o..o + nx];
                let sx_row = &sig_x[o..o + nx];
                let sy_row = &sig_y[o..o + nx];
                let prow = &mut prow[..nx];

                // i = 0 : vx[-1] = 0 (mur)
                let div0 = vxr[0] * inv_dx + (vyr[0] - vyp[0]) * inv_dy;
                let damping_0 = 1.0 / (1.0 + 0.5 * dt * (sx_row[0] + sy_row[0]));
                prow[0] = (prow[0] - dt * kr[0] * div0) * damping_0;

                for i in 1..nx {
                    let div = (vxr[i] - vxr[i - 1]) * inv_dx + (vyr[i] - vyp[i]) * inv_dy;
                    let damping = 1.0 / (1.0 + 0.5 * dt * (sx_row[i] + sy_row[i]));
                    prow[i] = (prow[i] - dt * kr[i] * div) * damping;
                }
            };
            for_each_row1(self.p.as_mut_slice(), nx, prs);
        }

        // 3. Enregistrement et avance du temps
        for (&(x, y), trace) in self.receivers.iter().zip(self.traces.iter_mut()) {
            trace.push(self.p.get(x, y));
        }
        self.time += dt as f64;
    }
}

// ---- Parcours des lignes : parallèle avec la feature `rayon`, séquentiel sinon.

#[cfg(feature = "rayon")]
fn for_each_row1<F>(a: &mut [f32], nx: usize, f: F)
where
    F: Fn(usize, &mut [f32]) + Sync + Send,
{
    use rayon::prelude::*;
    a.par_chunks_mut(nx).enumerate().for_each(|(y, r)| f(y, r));
}

#[cfg(not(feature = "rayon"))]
fn for_each_row1<F>(a: &mut [f32], nx: usize, f: F)
where
    F: Fn(usize, &mut [f32]) + Sync + Send,
{
    for (y, r) in a.chunks_mut(nx).enumerate() {
        f(y, r);
    }
}

#[cfg(feature = "rayon")]
fn for_each_row2<F>(a: &mut [f32], b: &mut [f32], nx: usize, f: F)
where
    F: Fn(usize, &mut [f32], &mut [f32]) + Sync + Send,
{
    use rayon::prelude::*;
    a.par_chunks_mut(nx)
        .zip(b.par_chunks_mut(nx))
        .enumerate()
        .for_each(|(y, (ra, rb))| f(y, ra, rb));
}

#[cfg(not(feature = "rayon"))]
fn for_each_row2<F>(a: &mut [f32], b: &mut [f32], nx: usize, f: F)
where
    F: Fn(usize, &mut [f32], &mut [f32]) + Sync + Send,
{
    for (y, (ra, rb)) in a.chunks_mut(nx).zip(b.chunks_mut(nx)).enumerate() {
        f(y, ra, rb);
    }
}
