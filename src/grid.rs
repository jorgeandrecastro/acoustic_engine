// SPDX-License-Identifier: GPL-2.0-or-later

/// Grille 2D de `f32` stockée en ligne-majeure (`index = y * nx + x`).
#[derive(Debug, Clone)]
pub struct Grid2D {
    pub nx: usize,
    pub ny: usize,
    pub dx: f32,
    pub dy: f32,
    data: Vec<f32>,
}

impl Grid2D {
    pub fn new(nx: usize, ny: usize, dx: f32, dy: f32) -> Self {
        Self {
            nx,
            ny,
            dx,
            dy,
            data: vec![0.0; nx * ny],
        }
    }

    #[inline(always)]
    pub fn index(&self, x: usize, y: usize) -> usize {
        y * self.nx + x
    }

    #[inline]
    pub fn get(&self, x: usize, y: usize) -> f32 {
        self.data[self.index(x, y)]
    }

    #[inline]
    pub fn set(&mut self, x: usize, y: usize, val: f32) {
        let idx = self.index(x, y);
        self.data[idx] = val;
    }

    pub fn fill(&mut self, val: f32) {
        self.data.fill(val);
    }

    /// Lecture seule des données brutes (export, rendu, boucles rapides).
    #[inline]
    pub fn as_slice(&self) -> &[f32] {
        &self.data
    }

    /// Écriture sur les données brutes.
    #[inline]
    pub fn as_mut_slice(&mut self) -> &mut [f32] {
        &mut self.data
    }

    /// Ligne `y` en lecture seule.
    #[inline]
    pub fn row(&self, y: usize) -> &[f32] {
        &self.data[y * self.nx..(y + 1) * self.nx]
    }
}
