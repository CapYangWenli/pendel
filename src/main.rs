const GRAV: f32 = 30.0;
const M1: f32 = 3.0;
const M2: f32 = 2.5;
const F_LENGTH: f32 = 2.5;
const S_LENGTH: f32 = 1.5;
const H: f32 = 0.01;
// const DEB: bool = false;

const L_COF: f32 = 100.0;


use macroquad::prelude::*;
// use std::{thread, time};

#[derive(Debug, Copy, Clone)]
struct Position {
    x: f32,
    y: f32
}

struct World {
    omega1: f32,
    omega2: f32,
    theta1: f32,
    theta2: f32,
}

impl World {
    fn total_energy(&self) -> f32 {
        let energy = 
        M1*(self.omega1*F_LENGTH).powi(2) / 2.0
        + M2*( 
                (self.omega2*S_LENGTH).powi(2) 
                + (self.omega1*F_LENGTH).powi(2)
                + 2.0*self.omega1*self.omega2*F_LENGTH*S_LENGTH*(self.theta1 - self.theta2).cos()
            ) / 2.0
        - M1*GRAV*F_LENGTH*self.theta1.cos()
        - M2*GRAV*(F_LENGTH*self.theta1.cos() + S_LENGTH*self.theta2.cos());

        energy
    }
    
}

#[macroquad::main("BasicSHapes")]
async fn main() {
    print!("Starting execution\n\n");


    let string_anchor: Position = Position { x: screen_width()/2.0, y:20.0};
    let mut world: World = World { omega1: 0.0, omega2: 0.0, theta1: 60.0, theta2: -10.0};
    let mut time = 0.0;

    let mut max_kinenergy: f32 = 0.0;

    world.theta1 = degtorad(&world.theta1);
    world.theta2 = degtorad(&world.theta2);

    
    loop {
        // print!("time: {}; (zeta1_old {} zeta2_old {} omega1_old {} omega2_old {}) -> ", time, world.theta1, world.theta2, world.omega1, world.omega2);
        update(&mut world, time);
        let energy = world.total_energy();
        if (max_kinenergy < energy) || (max_kinenergy == 0.0) {
            max_kinenergy = energy;
        }

        // print!("(zeta1_new {} zeta2_new {} omega1_new {} omega2_new {})\n", world.theta1, world.theta2, world.omega1, world.omega2);
        
        let fb_pos: Position = Position { x: string_anchor.x + F_LENGTH*L_COF*world.theta1.sin(), y: string_anchor.y + F_LENGTH*L_COF*world.theta1.cos() };
        let sb_pos: Position = Position { x: fb_pos.x + S_LENGTH*L_COF*world.theta2.sin(), y: fb_pos.y + S_LENGTH*L_COF*world.theta2.cos() };

        
        clear_background(BLACK);
        let _ball_pos_string= pos_string("Ball Position", &fb_pos);

        draw_line(string_anchor.x, string_anchor.y, fb_pos.x, fb_pos.y, 4.0, WHITE);
        draw_line(fb_pos.x, fb_pos.y, sb_pos.x, sb_pos.y, 4.0, WHITE);
        
        draw_circle(fb_pos.x, fb_pos.y, 15.0, WHITE);
        draw_circle(sb_pos.x, sb_pos.y, 15.0, WHITE);

        draw_text("IT WORKS!", 20.0, 20.0, 30.0, WHITE);
        draw_text(&format!("Max Total Energy: {}", max_kinenergy), screen_width() - 350.0, 20.0, 20.0, WHITE);
        draw_text(&format!("Cur Total Energy: {}", energy), screen_width() - 350.0, 40.0, 20.0, WHITE);


        time += H;
        // thread::sleep(time::Duration::from_millis(1));

        next_frame().await
    }
   
}

// tHeta = a
// omega = w
fn update(world: &mut World, time: f32) {

    //------- 1 ------
    // let thata1_n1 = rung_kutta(time, world.theta1, world, &der_zeta1);
    // let thata2_n1 = rung_kutta(time, world.theta2, world, &der_zeta2);
    
    // let omega1_n1 = rung_kutta(time, world.omega1, world, &der_omega1);
    // let omega2_n1 = rung_kutta(time, world.omega2, world, &der_omega2);    

    // (world.theta1, world.theta2, world.omega1, world.omega2) = 
    // (thata1_n1, thata2_n1, omega1_n1, omega2_n1)

    //------- 2 ------
    // let omega1_n1 = rung_kutta(time, world.omega1, world, &der_omega1);
    // let omega2_n1 = rung_kutta(time, world.omega2, world, &der_omega2); 
    // world.omega1 = omega1_n1;
    // world.omega2 = omega2_n1;

    // let thata1_n1 = rung_kutta(time, world.theta1, world, &der_zeta1);
    // let thata2_n1 = rung_kutta(time, world.theta2, world, &der_zeta2);

    // (world.theta1, world.theta2) = (thata1_n1, thata2_n1)

    //------- 3 ------
    world.omega1 = rung_kutta(time, world.omega1, world, &der_omega1);
    world.omega2 = rung_kutta(time, world.omega2, world, &der_omega2); 

    world.theta1 = rung_kutta(time, world.theta1, world, &der_zeta1);
    world.theta2 = rung_kutta(time, world.theta2, world, &der_zeta2);
    
}


fn rung_kutta(time: f32, var: f32, world: &World, fun: &dyn Fn(f32, f32, &World) -> f32) -> f32{
    
    let a: f32 = fun(time, var, world);
    // print!("a: {} is? var: {}\n", a, var);
    // b 1= f(t_n + H/2, w1_n + H*a/2)
    let b: f32 = fun(time + H/2.0, var + H*a/2.0, world);
    // c 1= f(t_n + H/2, w1_n + H*b/2)
    let c: f32 = fun(time + H/2.0, var + H*b/2.0, world);
    // d 1= f(t_n + H, w1_n + H*c)
    let d: f32 = fun(time + H, var + H*c, world);

    // let var = var + H*(a + 2.0*b + 2.0*c + d)/6.0;

    // let var = var + H*a/6.0 + H*b/3.0 + H*c/3.0 + H*d/6.0;

    
    // let new_var = var + cof;
    // print!("a: {} b: {} c:{} d:{} cof: {} var {} new_var {}\n", a, b, c, d, cof, var, new_var);
    // print!("Kutts cof: {} {}\n", cof, var);
    let var = var + fun(time + H, var, world)*H;

    var
}


fn der_omega1(_: f32, omega1: f32, world: &World) -> f32 {
    let theta1 = world.theta1;
    let theta2 = world.theta2;
    let omega2 = world.omega2;
    
    let delta_omega1: f32 = 
    (
        -GRAV*(2.0*M1 + M2)*theta1.sin() - 
        M2*GRAV*(theta1 - 2.0 * theta2).sin() -
        2.0*M2*omega2*omega2*S_LENGTH*(theta1 - theta2).sin() -
        M2*omega1*omega1*F_LENGTH*(2.0*(theta1 - theta2)).sin()
    ) /
    (F_LENGTH * ( 2.0*M1 + M2 - M2*(2.0*theta1 - 2.0*theta2).cos())); 

    delta_omega1
}

fn der_omega2(_: f32, omega2: f32, world: &World) -> f32 {
    let theta1 = world.theta1;
    let theta2 = world.theta2;
    let omega1 = world.omega1;
    
    let delta_omega2: f32 = 
    (
        2.0*(theta1-theta2).sin()*
        (
            (omega1*omega1*F_LENGTH*(M1 + M2)) + 
            GRAV*(M1 + M2)*theta1.cos() +
            omega2*omega2*S_LENGTH*M2*(theta1 - theta2).cos()
        ) 
    ) / 
    (S_LENGTH* 
    (  2.0*M1 + M2 - M2*(2.0*theta1 - 2.0*theta2).cos())
    ); 
    delta_omega2
}

fn der_zeta1(_: f32, _:f32, world: &World) -> f32 {
    // print!("der zets\n");
    world.omega1
}

fn der_zeta2(_: f32, _:f32, world: &World) -> f32 {
    // print!("der zets\n");
    world.omega2
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

fn degtorad(ang: &f32) -> f32 { 
    ang * std::f32::consts::PI / 180.0
}

//------------ DEBUG ----------

