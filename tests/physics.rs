// SPDX-License-Identifier: GPL-2.0-or-later

use acoustic_engine::{AcousticSolver2D, Waveform};

/// Maximum de |p| ; renvoie `f32::INFINITY` si une valeur est NaN ou infinie.
fn max_abs(s: &AcousticSolver2D) -> f32 {
    s.p.as_slice().iter().fold(0.0f32, |m, v| {
        if v.is_finite() {
            m.max(v.abs())
        } else {
            f32::INFINITY
        }
    })
}

/// Milieu normalisé : c = 1, rho = 1, dx = 1.
fn unit_solver(n: usize) -> AcousticSolver2D {
    AcousticSolver2D::new(n, n, 1.0, 1.0, 1.0, 1.0)
}

#[test]
fn cfl_formula() {
    let s = unit_solver(10);
    let expected = 1.0 / 2.0f32.sqrt();
    assert!((s.max_stable_dt() - expected).abs() < 1e-6);
}

#[test]
fn symmetry_of_the_wave() {
    let n = 101;
    let c = n / 2;
    let mut s = unit_solver(n);
    s.inject_source(c, c, 1.0);
    let dt = s.stable_dt(0.8);
    for _ in 0..60 {
        s.step(dt);
    }
    for d in [3usize, 10, 25] {
        let r = s.p.get(c + d, c);
        let l = s.p.get(c - d, c);
        let u = s.p.get(c, c + d);
        let b = s.p.get(c, c - d);
        let tol = 1e-4 * max_abs(&s);
        assert!((r - l).abs() < tol, "gauche/droite d={d}: {r} vs {l}");
        assert!((r - u).abs() < tol, "x/y d={d}: {r} vs {u}");
        assert!((r - b).abs() < tol, "haut/bas d={d}: {r} vs {b}");
    }
}

#[test]
fn causality_no_signal_faster_than_sound() {
    let n = 201;
    let c = n / 2;
    let mut s = unit_solver(n);
    s.add_source(c, c, Waveform::Ricker { freq: 0.05 }, 1.0);
    let rx = s.add_receiver(c + 60, c);
    let dt = s.stable_dt(0.5);
    // Le Ricker est centré à t = 20 ; le front ne peut pas arriver avant ~ 60/c.
    // Avant t = 40, rien ne doit être mesuré (l'onde a parcouru < 40 cellules depuis t=0).
    let steps = (40.0 / dt) as usize;
    for _ in 0..steps {
        s.step(dt);
    }
    let early = s.receiver(rx).iter().fold(0.0f32, |m, v| m.max(v.abs()));
    assert!(early < 1e-6, "signal trop précoce: {early}");

    // Plus tard, le signal arrive bien.
    for _ in 0..(120.0 / dt) as usize {
        s.step(dt);
    }
    let later = s.receiver(rx).iter().fold(0.0f32, |m, v| m.max(v.abs()));
    assert!(later > 1e-3, "le signal n'est jamais arrivé: {later}");
}

#[test]
fn stable_below_cfl() {
    let mut s = unit_solver(80);
    s.inject_source(40, 40, 1.0);
    let dt = s.stable_dt(0.95);
    for _ in 0..3000 {
        s.step(dt);
    }
    let m = max_abs(&s);
    assert!(m <= 1.0, "instabilité: {m}");
}

#[test]
fn unstable_above_cfl() {
    let mut s = unit_solver(80);
    s.inject_source(40, 40, 1.0);
    let dt = 1.3 * s.max_stable_dt();
    for _ in 0..600 {
        s.step(dt);
    }
    let m = max_abs(&s);
    assert!(m > 10.0, "devrait diverger, m = {m}");
}

#[test]
fn absorbing_boundary_removes_energy() {
    let run = |absorb: bool| {
        let mut s = unit_solver(120);
        if absorb {
            s.set_absorbing_boundary(20, 0.15);
        }
        s.inject_source(60, 60, 1.0);
        let dt = s.stable_dt(0.8);
        for _ in 0..(400.0 / dt) as usize {
            s.step(dt);
        }
        s.p.as_slice().iter().map(|v| v * v).sum::<f32>()
    };
    let rigid = run(false);
    let sponge = run(true);
    assert!(sponge < 0.05 * rigid, "rigide = {rigid}, éponge = {sponge}");
}

#[test]
fn heterogeneous_medium_changes_arrival_time() {
    // Bande rapide (c = 2) entre source et récepteur : arrivée plus tôt.
    let arrival =
        |fast: bool| {
            let mut s = unit_solver(201);
            if fast {
                s.set_medium_rect(0, 0, 201, 201, 1.0, 2.0);
            }
            s.add_source(100, 100, Waveform::Ricker { freq: 0.05 }, 1.0);
            let rx = s.add_receiver(150, 100);
            let dt = 0.4;
            for _ in 0..(150.0 / dt) as usize {
                s.step(dt);
            }
            let tr = s.receiver(rx);
            let (imax, _) = tr.iter().enumerate().fold((0, 0.0f32), |a, (i, v)| {
                if v.abs() > a.1 { (i, v.abs()) } else { a }
            });
            imax as f32 * dt
        };
    assert!(arrival(true) < arrival(false) - 10.0);
}
