// SPDX-License-Identifier: GPL-2.0-or-later
//! Enregistre un signal de pression en un point et l'écrit dans `receiver.csv`.

use acoustic_engine::{AcousticSolver2D, Waveform};
use std::{fs::File, io::Write};

fn main() -> std::io::Result<()> {
    // Air : rho = 1.2 kg/m³, c = 340 m/s, cellules de 0,5 m
    let mut s = AcousticSolver2D::new(300, 300, 0.5, 0.5, 1.2, 340.0);
    s.set_absorbing_boundary(25, 0.1);
    s.add_source(150, 150, Waveform::Ricker { freq: 100.0 }, 1.0);
    let rx = s.add_receiver(210, 150); // à 30 m de la source

    let dt = s.stable_dt(0.9);
    let steps = (0.25 / dt) as usize; // 250 ms
    for _ in 0..steps {
        s.step(dt);
    }

    let mut f = File::create("receiver.csv")?;
    writeln!(f, "t_s,pression")?;
    for (i, p) in s.receiver(rx).iter().enumerate() {
        writeln!(f, "{:e},{:e}", (i + 1) as f32 * dt, p)?;
    }
    println!("{steps} échantillons écrits dans receiver.csv (dt = {dt:e} s)");
    Ok(())
}
