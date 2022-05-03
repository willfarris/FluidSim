extern crate sdl2;


use sdl2::event::Event;
use sdl2::keyboard::Keycode;
use sdl2::pixels::Color;
use sdl2::rect::Rect;
use std::thread;
use std::time::Duration;

const WINDOW_WIDTH: u32 = 1024;
const GRID_SIZE: u32 = 8;
const NUM_SQUARES: usize = (WINDOW_WIDTH/GRID_SIZE) as usize;
const ARRAY_SIZE: usize = (NUM_SQUARES+2)*(NUM_SQUARES+2);

fn idx(field: &[f32; ARRAY_SIZE], x: usize, y: usize) -> f32 {
    field[x + (NUM_SQUARES+2)*y]
}

fn set_idx(field: &mut [f32; ARRAY_SIZE], val: f32, x: usize, y: usize) {
    field[x + (NUM_SQUARES+2)*y] = val;
}
    
fn set_bnd(b: u32, field: &mut [f32; ARRAY_SIZE]) {
    for i in 1..=NUM_SQUARES {
        set_idx(field, if b == 1 {-idx(field, 1, i)} else {idx(field, 1, i)}, 0, i);
        set_idx(field, if b == 1 {-idx(field, NUM_SQUARES, i)} else {idx(field, NUM_SQUARES, i)}, NUM_SQUARES+1, i);
        set_idx(field, if b == 2 {-idx(field, i, 1)} else {idx(field, i, 1)}, i, 0);
        set_idx(field, if b == 2 {-idx(field, i, NUM_SQUARES)} else {idx(field, i, NUM_SQUARES)}, i, NUM_SQUARES+1);
    }
    set_idx(field, 0.5*(idx(field, 1, 0)+idx(field, 0, 1)), 0, 0);

    set_idx(field, 0.5*(idx(field, 1, NUM_SQUARES+1)+idx(field, 0, NUM_SQUARES)), 0, NUM_SQUARES+1);
    
    set_idx(field, 0.5*(idx(field, NUM_SQUARES, 0)+idx(field, NUM_SQUARES+1, 1)), NUM_SQUARES+1, 0);
    
    set_idx(field, 0.5*(idx(field, NUM_SQUARES, NUM_SQUARES+1)+idx(field, NUM_SQUARES+1, NUM_SQUARES)), NUM_SQUARES+1, NUM_SQUARES+1);
}

fn diffuse(b: u32, mut x: [f32; ARRAY_SIZE], x0: [f32; ARRAY_SIZE], diff: f32, dt: f32) -> ([f32; ARRAY_SIZE], [f32; ARRAY_SIZE]) {
    let a = dt * diff * (NUM_SQUARES*NUM_SQUARES) as f32;

    for _ in 0..20 {
        for i in 1..=NUM_SQUARES {
            for j in 1..=NUM_SQUARES {
                let val = (idx(&x0, i, j) + a * (idx(&x, i-1, j) + idx(&x, i+1, j) + idx(&x, i, j - 1) + idx(&x, i, j + 1))) / (1. + 4.*a);
                set_idx(&mut x, val, i, j);
            }
        }
        set_bnd(b, &mut x);
    }

    (x, x0)
}

fn advect(b: u32, mut d: [f32; ARRAY_SIZE], d0: [f32; ARRAY_SIZE], u: &[f32; ARRAY_SIZE], v: &[f32; ARRAY_SIZE], dt: f32) -> ([f32; ARRAY_SIZE], [f32; ARRAY_SIZE]) {

    let dt0 = dt*NUM_SQUARES as f32;
    for i in 1..NUM_SQUARES {
        for j in 1..NUM_SQUARES {
            let mut x = i as f32 - dt0*idx(&u, i, j);
            let mut y = j as f32 - dt0*idx(&v, i, j);

            if x < 0.5 {
                x = 0.5;
            }
            if x > NUM_SQUARES as f32 + 0.5 {
                x = NUM_SQUARES as f32 + 0.5;
            }
            let i0 = x.floor() as usize;
            let i1 = i0 + 1;

            if y < 0.5 {
                y = 0.5;
            }
            if y > NUM_SQUARES as f32 + 0.5 {
                y = NUM_SQUARES as f32 + 0.5;
            }
            let j0 = y.floor() as usize;
            let j1 = j0 + 1;

            let s1 = x - i0 as f32;
            let s0 = 1. - s1;

            let t1 = y - j0 as f32;
            let t0 = 1. - t1;

            let d_new = s0*(t0*idx(&d0, i0, j0) + t1*idx(&d0, i0, j1)) + s1*(t0*idx(&d0, i1, j0) + t1*idx(&d0, i1, j1));
            set_idx(&mut d, d_new, i, j);
        }
    }
    set_bnd(b, &mut d);
    (d, d0)
}

fn project(mut u: [f32; ARRAY_SIZE], mut v: [f32; ARRAY_SIZE]) -> ([f32; ARRAY_SIZE], [f32; ARRAY_SIZE]) {
    let h = 1./NUM_SQUARES as f32;
    let mut div = [0f32; ARRAY_SIZE];
    let mut p = [0f32; ARRAY_SIZE];

    for i in 1..=NUM_SQUARES {
        for j in 1..=NUM_SQUARES {
            let div_val = -0.5 * h * (idx(&u, i+1, j) - idx(&u, i-1, j) + idx(&v, i, j+1) - idx(&v, i, j-1));
            set_idx(&mut div, div_val, i, j);
        }
    }
    set_bnd(0, &mut div);
    set_bnd(0, &mut p);

    for _ in 0..20 {
        for i in 1..NUM_SQUARES {
            for j in 1..NUM_SQUARES {
                let p_val = (idx(&div, i, j) + idx(&p, i-1, j) + idx(&p, i+1, j) + idx(&p, i, j-1) + idx(&p, i, j+1))/4.;
                set_idx(&mut p, p_val, i, j);
            }
        }
        set_bnd(0, &mut p);
    }

    for i in 1..NUM_SQUARES {
        for j in 1..NUM_SQUARES {
            let u_val = idx(&u, i, j) - 0.5 * (idx(&p, i+1, j) - idx(&p, i-1, j))/h;
            set_idx(&mut u, u_val, i, j);

            let v_val = idx(&v, i, j) - 0.5 * (idx(&p, 1, j+1) - idx(&p, i, j-1))/h;
            set_idx(&mut v, v_val, i, j);
        }
    }
    set_bnd(1, &mut u);
    set_bnd(2, &mut v);

    (u, v)
}

fn dens_step(mut x: [f32; ARRAY_SIZE], mut x0: [f32; ARRAY_SIZE], mut u: [f32; ARRAY_SIZE], mut v: [f32; ARRAY_SIZE], diff: f32, dt: f32) -> ([f32; ARRAY_SIZE], [f32; ARRAY_SIZE], [f32; ARRAY_SIZE], [f32; ARRAY_SIZE]) {
    
    let temp = x;
    x = x0;
    x0 = temp;

    (x, x0) = diffuse(0, x, x0, diff, dt);
    
    let temp = x;
    x = x0;
    x0 = temp;
    
    (x, x0) = advect(0, x, x0, &u, &v, dt);
    
    (x, x0, u, v)
}

fn vel_step(mut u: [f32; ARRAY_SIZE], mut v: [f32; ARRAY_SIZE], mut u0: [f32; ARRAY_SIZE], mut v0: [f32; ARRAY_SIZE], visc: f32, dt: f32 ) -> ([f32; ARRAY_SIZE], [f32; ARRAY_SIZE], [f32; ARRAY_SIZE], [f32; ARRAY_SIZE]) {
    let temp = u;
    u = u0;
    u0 = temp;
    (u, u0) = diffuse(1, u, u0, visc, dt);

    let temp = v;
    v = v0;
    v0 = temp;
    (v, v0) = diffuse(2, v, v0, visc, dt);

    (u, v) = project(u, v);

    let temp = u;
    u = u0;
    u0 = temp;

    let temp = v;
    v = v0;
    v0 = temp;

    (u, u0) = advect(1, u, u0, &u0, &v0, 0.01);
    (v, v0) = advect(2, v, v0, &u0, &v0, 0.01);

    (u, v) = project(u, v);

    (u, v, u0, v0)
}

fn main() {
    let sdl_context = sdl2::init().unwrap();
    let video_subsystem = sdl_context.video().unwrap();
    let window = video_subsystem
        .window("Fluid Simulation", WINDOW_WIDTH as u32, WINDOW_WIDTH as u32)
        .position_centered()
        .build()
        .map_err(|e| e.to_string())
        .unwrap();
    println!("Hello, world!");

    let mut canvas = window
        .into_canvas()
        .software()
        .build()
        .map_err(|e| e.to_string()).unwrap();
        
    
    let mut u_fluid: [f32; ARRAY_SIZE] = [0f32; ARRAY_SIZE];
    let mut v_fluid: [f32; ARRAY_SIZE] = [0f32; ARRAY_SIZE];
    let mut d_fluid: [f32; ARRAY_SIZE] = [0f32; ARRAY_SIZE];

    let mut u_prev_fluid: [f32; ARRAY_SIZE] = [0f32; ARRAY_SIZE];
    let mut v_prev_fluid: [f32; ARRAY_SIZE] = [0f32; ARRAY_SIZE];
    let mut d_prev_fluid: [f32; ARRAY_SIZE] = [0f32; ARRAY_SIZE];

    /*
    for i in 1..NUM_SQUARES {
        set_idx(&mut u_fluid, 20.* i as f32 / NUM_SQUARES as f32, i, NUM_SQUARES/2);
    }
    */

    let mut adding = false;
    let mut add_x = 0;
    let mut add_y = 0;

    'mainloop: loop {
        thread::sleep(Duration::from_millis(20));
        for event in sdl_context.event_pump().unwrap().poll_iter() {
            match event {
                Event::KeyDown {
                    keycode: Some(Keycode::Escape),
                    ..
                } | Event::Quit {..} => break 'mainloop,
                Event::MouseButtonDown {x, y, ..} => {
                    adding = !adding;
                    add_x = x;
                    add_y = y;
                },
                _ => {}
            }
        }

        if adding {
            set_idx(&mut d_fluid, 10., (add_x as u32 / GRID_SIZE) as usize, (add_y as u32 / GRID_SIZE) as usize);
            //set_idx(&mut u_fluid, 100., (x as u32 / GRID_SIZE) as usize, (y as u32 / GRID_SIZE) as usize);
            set_idx(&mut v_fluid, -100., (add_x as u32 / GRID_SIZE) as usize, (add_y as u32 / GRID_SIZE) as usize);
        }

        let visc = 0.0001;
        let dt = 1./144.;

        // Velocity step
        (u_fluid, v_fluid, u_prev_fluid, v_prev_fluid) = vel_step(u_fluid, v_fluid, u_prev_fluid, v_prev_fluid, visc, dt);

        // Density step
        (d_fluid, d_prev_fluid, u_fluid, v_fluid) = dens_step(d_fluid, d_prev_fluid, u_fluid, v_fluid, visc, dt);
        

        
        canvas.set_draw_color(Color::RGBA(0, 0, 0, 255));
        canvas.clear();
        for i in 1..NUM_SQUARES+2 {
            for j in 1..NUM_SQUARES+2 {
                let red = 0;// (idx(&mut &u_fluid, i, j) * 255.0).abs().floor() as u8;
                let green = 0;// (idx(&mut v_fluid, i, j) * 255.0).abs().floor() as u8;
                let blue = (idx(&d_fluid, i, j) * 512.0).floor() as u8;
                canvas.set_draw_color(Color::RGBA(red, green, blue, 255));
                canvas
                    .fill_rect(Rect::new((i as u32 * GRID_SIZE) as i32, (j as u32 * GRID_SIZE) as i32, GRID_SIZE, GRID_SIZE))
                    .expect("could not fill rect");
                
                /*
                let start_x = (i as u32 * GRID_SIZE + GRID_SIZE/2) as i32;
                let start_y = (j as u32 * GRID_SIZE + GRID_SIZE/2) as i32;
                let end_x = start_x + (idx(&mut &u_fluid, i, j) * GRID_SIZE as f32).floor() as i32;
                let end_y = start_y + (idx(&mut v_fluid, i, j) * GRID_SIZE as f32).floor() as i32;
                canvas.set_draw_color(Color::RGBA(122, 122, 122, 255));
                canvas.draw_line((start_x, start_y), (end_x, end_y)).expect("Could not draw arrow");
                */
            }
        }
        
        canvas.present();
    }
}
