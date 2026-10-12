// SPDX-License-Identifier: GPL-2.0-or-later
//! Mesure le débit du solveur : cargo run --release --example bench [--features rayon]

use acoustic_engine::AcousticSolver2D;
use std::time::Instant;

fn main() {
    let n = 1024;
    let mut s = AcousticSolver2D::new(n, n, 1.0, 1.0, 1.0, 1.0);
    s.inject_source(n / 2, n / 2, 1.0);
    let dt = s.stable_dt(0.9);

    for _ in 0..10 {
        s.step(dt); // échauffement
    }
    let steps = 200;
    let t0 = Instant::now();
    for _ in 0..steps {
        s.step(dt);
    }
    let secs = t0.elapsed().as_secs_f64();
    let mcells = (n * n * steps) as f64 / secs / 1e6;
    println!("{n}x{n}, {steps} pas : {secs:.2} s  ->  {mcells:.0} Mcellules/s");
}
