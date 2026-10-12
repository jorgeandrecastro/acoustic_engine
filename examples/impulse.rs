// SPDX-License-Identifier: GPL-2.0-or-later

use acoustic_engine::AcousticSolver2D;

fn main() {
    // Grille 100x100, pas 1 m, eau : rho = 1000 kg/m³, c = 1500 m/s
    let mut solver = AcousticSolver2D::new(100, 100, 1.0, 1.0, 1000.0, 1500.0);
    solver.set_absorbing_boundary(15, 0.1);

    let dt = solver.stable_dt(0.9);
    println!("dt = {dt:e} s (limite CFL respectée)");

    solver.inject_source(50, 50, 1000.0);
    for _ in 0..100 {
        solver.step(dt);
    }

    println!("t = {:.3e} s, pression au centre : {}", solver.time(), solver.p.get(50, 50));
}
