// SPDX-License-Identifier: GPL-2.0-or-later
//! Visualisation temps réel : cargo run --release --example viewer [--features rayon]
//! Rouge = surpression, bleu = dépression, gris = obstacle. Échap pour quitter.

use acoustic_engine::{AcousticSolver2D, Waveform};
use minifb::{Key, Window, WindowOptions};
use std::time::Duration;

const N: usize = 300;
const SCALE: usize = 3;
const STEPS_PER_FRAME: usize = 4;

fn main() {
    // Air : rho = 1.2 kg/m³, c = 340 m/s, cellules de 1 m
    let mut s = AcousticSolver2D::new(N, N, 1.0, 1.0, 1.2, 340.0);
    s.set_absorbing_boundary(25, 0.1);

    // Obstacle : mur très dense (même c, pour ne pas réduire le pas de temps)
    s.set_medium_rect(190, 90, 200, 210, 1200.0, 340.0);

    s.add_source(100, 150, Waveform::Ricker { freq: 25.0 }, 1.0);
    let dt = s.stable_dt(0.9);

    let (w, h) = (N * SCALE, N * SCALE);
    let mut buffer = vec![0u32; w * h];
    let mut window = Window::new("acoustic-engine", w, h, WindowOptions::default())
        .expect("impossible de créer la fenêtre");
    window.limit_update_rate(Some(Duration::from_micros(16_600)));

    let mut gain = 1e-6f32;
    while window.is_open() && !window.is_key_down(Key::Escape) {
        for _ in 0..STEPS_PER_FRAME {
            s.step(dt);
        }

        // Gain automatique (décroît lentement, remonte avec le pic courant)
        let peak = s.p.as_slice().iter().fold(0.0f32, |m, v| m.max(v.abs()));
        gain = (gain * 0.995).max(peak).max(1e-6);

        for y in 0..N {
            for x in 0..N {
                let color = if s.rho().get(x, y) > 100.0 {
                    0x00404040
                } else {
                    let t = (s.p.get(x, y) / gain).clamp(-1.0, 1.0);
                    let fade = ((1.0 - t.abs()) * 255.0) as u32;
                    if t >= 0.0 {
                        0x00FF0000 | (fade << 8) | fade
                    } else {
                        (fade << 16) | (fade << 8) | 0x000000FF
                    }
                };
                for dy in 0..SCALE {
                    let row = (y * SCALE + dy) * w + x * SCALE;
                    buffer[row..row + SCALE].fill(color);
                }
            }
        }
        window.update_with_buffer(&buffer, w, h).unwrap();
    }
}
