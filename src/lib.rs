// Copyright (C) 2026 Ton Nom
// SPDX-License-Identifier: GPL-2.0-or-later

//! Solveur acoustique 2D par différences finies (FDTD) sur grille décalée.
//!
//! ```
//! use acoustic_engine::{AcousticSolver2D, Waveform};
//!
//! let mut s = AcousticSolver2D::new(64, 64, 1.0, 1.0, 1000.0, 1500.0);
//! s.set_absorbing_boundary(8, 0.1);
//! s.add_source(32, 32, Waveform::Ricker { freq: 100.0 }, 1.0);
//! let rx = s.add_receiver(10, 32);
//! let dt = s.stable_dt(0.9);
//! for _ in 0..50 {
//!     s.step(dt);
//! }
//! assert_eq!(s.receiver(rx).len(), 50);
//! ```

pub mod grid;
pub mod solver;

pub use grid::Grid2D;
pub use solver::{AcousticSolver2D, Waveform};
