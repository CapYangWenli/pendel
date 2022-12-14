const GRAV: f32 = 30.0;

const M1: f32 = 4.0;
const M2: f32 = 2.5;
const L1: f32 = 2.0;
const L2: f32 = 1.5;
const H: f32 = 0.01;
// const DEB: bool = false;

const L_COF: f32 = 100.0;
const THETA_1: f32 = 60.0;
const THETA_2: f32 = -30.0;

const OMEGA_1: f32 = 0.0;
const OMEGA_2: f32 = 0.0;

use macroquad::{miniquad::date::now, prelude::*};

#[derive(Debug, Copy, Clone)]
struct Position {
    x: f32,
    y: f32,
}

fn total_energy(var_vec: &Vec<f32>) -> f32 {
    let theta1 = var_vec[0];
    let theta2 = var_vec[1];
    let omega1 = var_vec[2];
    let omega2 = var_vec[3];

    let energy = M1 * (omega1 * L1).powi(2) / 2.0
        + M2 * ((omega2 * L2).powi(2)
            + (omega1 * L1).powi(2)
            + 2.0 * omega1 * omega2 * L1 * L2 * (theta1 - theta2).cos())
            / 2.0
        - M1 * GRAV * L1 * theta1.cos()
        - M2 * GRAV * (L1 * theta1.cos() + L2 * theta2.cos());

    energy
}

#[macroquad::main("Simple Pendel")]
async fn main() {
    print!("Starting execution\n\n");

    // var_vec[0] = theta_1
    // var_vec[1] = theta_2
    // var_vec[2] = omega_1
    // var_vec[3] = omega_2
    // var_vec[4] = time
    let mut var_vec: Vec<f32> = vec![degtorad(THETA_1), degtorad(THETA_2), OMEGA_1, OMEGA_2, 0.0];
    let fn_vec: Vec<&dyn Fn(&Vec<f32>) -> f32> =
        vec![&der_zeta1, &der_zeta2, &der_omega1, &der_omega2, &der_time];

    let mut max_kinenergy: f32 = 0.0;
    let mut traj: Vec<(Position, Position)> = Vec::new();
    let mut frames: i32 = 0;

    loop {
        let string_anchor: Position = Position {
            x: screen_width() / 2.0,
            y: 20.0,
        };
        update(&mut var_vec, &fn_vec);
        let theta_1 = &var_vec[0];
        let theta_2 = &var_vec[1];

        let energy = total_energy(&var_vec);
        if (max_kinenergy < energy) || (max_kinenergy == 0.0) {
            max_kinenergy = energy;
        }


        let fb_pos: Position = Position {
            x: string_anchor.x + L1 * L_COF * theta_1.sin(),
            y: string_anchor.y + L1 * L_COF * theta_1.cos(),
        };
        let sb_pos: Position = Position {
            x: fb_pos.x + L2 * L_COF * theta_2.sin(),
            y: fb_pos.y + L2 * L_COF * theta_2.cos(),
        };

        clear_background(BLACK);
        let _ball_pos_string = pos_string("Ball Position", &fb_pos);

        draw_text("IT WORKS!", 20.0, 20.0, 30.0, WHITE);
        draw_text(
            &format!("FPS: {}", get_fps()),
            screen_width() - 350.0,
            20.0,
            20.0,
            WHITE,
        );
        draw_text(
            &format!("Cur Total Energy: {}", energy),
            screen_width() - 350.0,
            40.0,
            20.0,
            WHITE,
        );

        traj.push((fb_pos.clone(), sb_pos.clone()));

        if traj.len() > 60 * 3 {
            traj.remove(0);
        }

        let mut last_f_ball: Position = Position { x: 0.0, y: 0.0 };
        let mut last_s_ball: Position = Position { x: 0.0, y: 0.0 };

        for (i, (first_ball, sec_ball)) in traj.iter().enumerate() {
            draw_line(
                first_ball.x,
                first_ball.y,
                last_f_ball.x,
                last_f_ball.y,
                4.0 * (i as f32 / traj.len() as f32),
                BLUE,
            );
            draw_line(
                sec_ball.x,
                sec_ball.y,
                last_s_ball.x,
                last_s_ball.y,
                4.0 * (i as f32 / traj.len() as f32),
                RED,
            );
            last_f_ball = first_ball.clone();
            last_s_ball = sec_ball.clone();
        }

        draw_line(
            string_anchor.x,
            string_anchor.y,
            fb_pos.x,
            fb_pos.y,
            4.0,
            WHITE,
        );
        draw_line(fb_pos.x, fb_pos.y, sb_pos.x, sb_pos.y, 4.0, WHITE);

        draw_circle(fb_pos.x, fb_pos.y, 15.0, BLUE);
        draw_circle(sb_pos.x, sb_pos.y, 15.0, RED);

        if frames > 60 * 100 {
            frames = 0;
        }
        frames += 1;

        next_frame().await
    }
}

fn update(var_vec: &mut Vec<f32>, fn_vec: &Vec<&dyn Fn(&Vec<f32>) -> f32>) {
    let new_var_vec = rung_kutta(var_vec, fn_vec);
    *var_vec = new_var_vec;
}

fn rung_kutta(var_vec: &Vec<f32>, fn_vec: &Vec<&dyn Fn(&Vec<f32>) -> f32>) -> Vec<f32> {
    // var_vec[0] = theta_1
    // var_vec[1] = theta_2
    // var_vec[2] = omega_1
    // var_vec[3] = omega_2
    // var_vec[4] = time

    let a: Vec<f32> = comp(var_vec, fn_vec);

    // b 1= f(t_n + H/2, w1_n + H*a/2)

    let a_h = (&a).into_iter().map(|x| H / 2.0 * x).collect::<Vec<_>>();
    let b: Vec<f32> = comp(&add(&a_h, var_vec), fn_vec);

    // c 1= f(t_n + H/2, w1_n + H*b/2)
    let b_h = (&b).into_iter().map(|x| H / 2.0 * x).collect::<Vec<_>>();
    let c: Vec<f32> = comp(&add(&b_h, &var_vec), fn_vec);

    // d 1= f(t_n + H, w1_n + H*c)
    let c_h = (&c).into_iter().map(|x| H * x).collect::<Vec<_>>();
    let d: Vec<f32> = comp(&add(&c_h, &var_vec), fn_vec);

    //  var = var + H*(a + 2.0*b + 2.0*c + d)/6.0;
    let b_2 = (&b).into_iter().map(|x| 2.0 * x).collect::<Vec<_>>();
    let c_2 = (&c).into_iter().map(|x| 2.0 * x).collect::<Vec<_>>();

    let coff = add(&add(&add(&a, &b_2), &c_2), &d);
    let coff_h = coff.into_iter().map(|x| H * x / 6.0).collect::<Vec<_>>();

    let var = add(&var_vec, &coff_h);

    var
}

fn comp(var_vec: &Vec<f32>, fn_vec: &Vec<&dyn Fn(&Vec<f32>) -> f32>) -> Vec<f32> {
    let mut ans: Vec<f32> = Vec::new();

    for fnc in fn_vec {
        ans.push(fnc(var_vec));
    }

    ans
}

fn der_omega1(world: &Vec<f32>) -> f32 {
    let theta1 = world[0];
    let theta2 = world[1];
    let omega1 = world[2];
    let omega2 = world[3];

    let delta_omega1: f32 = (-GRAV * (2.0 * M1 + M2) * theta1.sin()
        - M2 * GRAV * (theta1 - 2.0 * theta2).sin()
        - 2.0 * M2 * omega2 * omega2 * L2 * (theta1 - theta2).sin()
        - M2 * omega1 * omega1 * L1 * (2.0 * (theta1 - theta2)).sin())
        / (L1 * (2.0 * M1 + M2 - M2 * (2.0 * theta1 - 2.0 * theta2).cos()));

    delta_omega1
}

fn der_omega2(world: &Vec<f32>) -> f32 {
    let theta1 = world[0];
    let theta2 = world[1];
    let omega1 = world[2];
    let omega2 = world[3];

    let delta_omega2: f32 = (2.0
        * (theta1 - theta2).sin()
        * ((omega1 * omega1 * L1 * (M1 + M2))
            + GRAV * (M1 + M2) * theta1.cos()
            + omega2 * omega2 * L2 * M2 * (theta1 - theta2).cos()))
        / (L2 * (2.0 * M1 + M2 - M2 * (2.0 * theta1 - 2.0 * theta2).cos()));
    delta_omega2
}

fn der_zeta1(world: &Vec<f32>) -> f32 {
    let omega1 = world[2];

    omega1
}

fn der_zeta2(world: &Vec<f32>) -> f32 {
    let omega2 = world[3];

    omega2
}

fn der_time(_: &Vec<f32>) -> f32 {
    1.0
}

fn add(a: &Vec<f32>, b: &Vec<f32>) -> Vec<f32> {
    let z = a
        .iter()
        .zip(b.iter())
        .map(|(&b, &v)| b + v)
        .collect::<Vec<_>>();
    z
}

fn pos_string(descr: &str, pos: &Position) -> String {
    let mut pos_string = String::from(descr);
    pos_string.push_str(" (x: ");
    pos_string.push_str(&pos.x.to_string());
    pos_string.push_str(" y: ");
    pos_string.push_str(&pos.y.to_string());
    pos_string.push_str(")");

    pos_string
}

fn degtorad(ang: f32) -> f32 {
    ang * std::f32::consts::PI / 180.0
}

