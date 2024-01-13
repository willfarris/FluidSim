extern crate sdl2;


use sdl2::event::Event;
use sdl2::keyboard::Keycode;
use sdl2::pixels::Color;
use sdl2::rect::Rect;

const WINDOW_WIDTH: u32 = 1024;
const GRID_SIZE: u32 = 4;
const NUM_SQUARES: usize = (WINDOW_WIDTH/GRID_SIZE) as usize;
const ARRAY_SIZE: usize = (NUM_SQUARES+2)*(NUM_SQUARES+2);

struct Fluid {
    visc: f32,

    #[allow(dead_code)]
    s: [f32; ARRAY_SIZE],

    d: [f32; ARRAY_SIZE],
    u: [f32; ARRAY_SIZE],
    v: [f32; ARRAY_SIZE],

    d0: [f32; ARRAY_SIZE],
    u0: [f32; ARRAY_SIZE],
    v0: [f32; ARRAY_SIZE],

    _div: [f32; ARRAY_SIZE],
    _p: [f32; ARRAY_SIZE],
}

trait FluidArray<T> {
    fn idx(&self, x: usize, y: usize) -> T;
    fn set_idx(&mut self, val: T, x: usize, y: usize);
    fn set_bnd(&mut self, b: u32);
    fn diffuse(&mut self, b: u32, x0: &[T; ARRAY_SIZE], diff: T, dt: f32);
    fn advect(&mut self, b: u32, d0: &[T; ARRAY_SIZE], u: &[T; ARRAY_SIZE], v: &[T; ARRAY_SIZE], dt: f32);
    //fn project(u: &mut [T; ARRAY_SIZE], v: &mut [T; ARRAY_SIZE], p: &mut [T; ARRAY_SIZE], div: &mut [T; ARRAY_SIZE]) -> ([f32; ARRAY_SIZE], [f32; ARRAY_SIZE]);
}

impl FluidArray<f32> for [f32; ARRAY_SIZE] {
    fn idx(&self, x: usize, y: usize) -> f32 {
        self[x + (NUM_SQUARES+2)*y]
    }

    fn set_idx(&mut self, val: f32, x: usize, y: usize) {
        self[x + (NUM_SQUARES+2)*y] = val;
    }

    fn set_bnd(&mut self, b: u32) {
        for i in 1..=NUM_SQUARES {
            self.set_idx( if b == 1 {-self.idx( 1, i)} else {self.idx( 1, i)}, 0, i);
            self.set_idx(if b == 1 {-self.idx(NUM_SQUARES, i)} else {self.idx(NUM_SQUARES, i)}, NUM_SQUARES+1, i);
            self.set_idx(if b == 2 {-self.idx(i, 1)} else {self.idx(i, 1)}, i, 0);
            self.set_idx(if b == 2 {-self.idx(i, NUM_SQUARES)} else {self.idx(i, NUM_SQUARES)}, i, NUM_SQUARES+1);
        }
        self.set_idx(0.5*(self.idx(1, 0)+self.idx(0, 1)), 0, 0);
        self.set_idx(0.5*(self.idx(1, NUM_SQUARES+1)+self.idx(0, NUM_SQUARES)), 0, NUM_SQUARES+1);
        self.set_idx(0.5*(self.idx(NUM_SQUARES, 0)+self.idx(NUM_SQUARES+1, 1)), NUM_SQUARES+1, 0);
        self.set_idx(0.5*(self.idx(NUM_SQUARES, NUM_SQUARES+1)+self.idx(NUM_SQUARES+1, NUM_SQUARES)), NUM_SQUARES+1, NUM_SQUARES+1);
    }

    fn diffuse(&mut self, b: u32, x0: &[f32; ARRAY_SIZE], diff: f32, dt: f32) {
        let a = dt * diff * (NUM_SQUARES*NUM_SQUARES) as f32;
    
        for _ in 0..20 {
            for i in 1..=NUM_SQUARES {
                for j in 1..=NUM_SQUARES {
                    let val = (x0.idx(i, j) + a * (self.idx(i-1, j) + self.idx(i+1, j) + self.idx(i, j - 1) + self.idx(i, j + 1))) / (1. + 4.*a);
                    self.set_idx( val, i, j);
                }
            }
            self.set_bnd(b,);
        }
    }

    fn advect(&mut self, b: u32, d0: &[f32; ARRAY_SIZE], u: &[f32; ARRAY_SIZE], v: &[f32; ARRAY_SIZE], dt: f32) {

        let dt0 = dt*NUM_SQUARES as f32;

        for i in 1..NUM_SQUARES {
            for j in 1..NUM_SQUARES {
                let mut x = i as f32 - dt0*u.idx(i, j);
                let mut y = j as f32 - dt0*v.idx(i, j);
    
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
    
                let d_new = s0*(t0*d0.idx(i0, j0) + t1*d0.idx(i0, j1)) + s1*(t0*d0.idx(i1, j0) + t1*d0.idx(i1, j1));
                self.set_idx(d_new, i, j);
            }
        }
        self.set_bnd(b);
    }
}

impl Fluid {
    fn new(visc: f32) -> Self {
        Self {
            visc,

            s: [0f32; ARRAY_SIZE],

            d: [0f32; ARRAY_SIZE],
            u: [0f32; ARRAY_SIZE],
            v: [0f32; ARRAY_SIZE],

            d0: [0f32; ARRAY_SIZE],
            u0: [0f32; ARRAY_SIZE],
            v0: [0f32; ARRAY_SIZE],

            _div: [0f32; ARRAY_SIZE],
            _p: [0f32; ARRAY_SIZE],
        }
    }

    fn project(&mut self) {//u: &mut [f32; ARRAY_SIZE], v: &mut [f32; ARRAY_SIZE], p: &mut [f32; ARRAY_SIZE], div: &mut [f32; ARRAY_SIZE]) -> ([f32; ARRAY_SIZE], [f32; ARRAY_SIZE]) {
        let h = 1./NUM_SQUARES as f32;
        //let mut div = [0f32; ARRAY_SIZE];
        //let mut p = [0f32; ARRAY_SIZE];
    
        for i in 1..=NUM_SQUARES {
            for j in 1..=NUM_SQUARES {
                let div_val = -0.5 * h * (self.u.idx(i+1, j) - self.u.idx(i-1, j) + self.v.idx(i, j+1) - self.v.idx(i, j-1));
                self._div.set_idx(div_val, i, j);
                self._p.set_idx(0f32, i, j);
            }
        }
        self._div.set_bnd(0);
        self._p.set_bnd(0);
    
        for _ in 0..20 {
            for i in 1..NUM_SQUARES {
                for j in 1..NUM_SQUARES {
                    let p_val = (self._div.idx(i, j) + self._p.idx(i-1, j) + self._p.idx(i+1, j) + self._p.idx( i, j-1) + self._p.idx( i, j+1))/4.;
                    self._p.set_idx(p_val, i, j);
                }
            }
            self._p.set_bnd(0);
        }
    
        for i in 1..NUM_SQUARES {
            for j in 1..NUM_SQUARES {
                let u_val = self.u.idx(i, j) - 0.5 * (self._p.idx( i+1, j) - self._p.idx(i-1, j))/h;
                self.u.set_idx(u_val, i, j);
    
                let v_val = self.v.idx(i, j) - 0.5 * (self._p.idx(1, j+1) - self._p.idx( i, j-1))/h;
                self.v.set_idx(v_val, i, j);
            }
        }
        self.u.set_bnd(1);
        self.v.set_bnd(2);
    }

    fn dens_step(&mut self, dt: f32) {
        std::mem::swap(&mut self.d, &mut self.d0);
        self.d.diffuse(0, &self.d0, self.visc, dt);
        std::mem::swap(&mut self.d, &mut self.d0);
        self.d.advect(0, &self.d0, &self.u, &self.v, dt);
    }
    
    fn vel_step(&mut self, dt: f32) {
        std::mem::swap(&mut self.u, &mut self.u0);
        self.u.diffuse(1, &self.u0, self.visc, dt);
        std::mem::swap(&mut self.v, &mut self.v0);
        self.v.diffuse(2, &self.v0, self.visc, dt);
        self.project();
        std::mem::swap(&mut self.u, &mut self.u0);
        std::mem::swap(&mut self.v, &mut self.v0);
        self.u.advect(1, &self.u0, &self.u0, &self.v0, 0.01);
        self.v.advect(2, &self.v0, &self.u0, &self.v0, 0.01);
        self.project();
    }

    fn step_simulation(&mut self, adding_dense: bool, dt: f32) {
        let density_added = 20.;
        
        // Add fluid
        if adding_dense {
            self.d.set_idx(density_added, NUM_SQUARES / 2, NUM_SQUARES / 2);
            self.d0.set_idx(density_added, NUM_SQUARES / 2, NUM_SQUARES / 2);
        }

        self.vel_step(dt);
        self.dens_step(dt);
    }
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

    let mut canvas = window
        .into_canvas()
        .software()
        .build()
        .map_err(|e| e.to_string()).unwrap();

    let visc = 0.005;
    let dt = 1./144.;

    let mut fluid = Box::new(Fluid::new(visc));

    let mut adding_vel = false;
    let mut adding_dens = false;
    let mut show_vec_field = false;

    'mainloop: loop {
        
        for event in sdl_context.event_pump().unwrap().poll_iter() {
            match event {
                Event::KeyDown {
                    keycode: Some(Keycode::Escape),
                    ..
                } | Event::Quit {..} => break 'mainloop,
                Event::KeyDown {
                    keycode: Some(Keycode::Space),
                    ..
                } => { adding_dens = !adding_dens; },
                Event::KeyDown {
                    keycode: Some(Keycode::V),
                    ..
                } => { show_vec_field = !show_vec_field; },
                Event::MouseButtonDown {..} => {
                    adding_vel = true;
                },
                Event::MouseButtonUp {..} => {
                    adding_vel = false;
                },
                Event::MouseMotion { x, y, xrel, yrel, ..} => {
                    let scale = 10.0;
                    let x = x as usize / GRID_SIZE as usize;
                    let y = y as usize / GRID_SIZE as usize;
                    
                    if adding_vel {
                        fluid.u.set_idx(scale * xrel as f32, x, y);
                        fluid.v.set_idx(scale * yrel as f32, x, y);

                        fluid.u0.set_idx(scale * xrel as f32, x , y);
                        fluid.v0.set_idx(scale * yrel as f32, x, y);
                        for i in 0..4 {
                            for j in 0..4 {
                                let x = (x+i);
                                let y = (y+j);
                                fluid.d.set_idx(1.0, x, y);
                                fluid.d0.set_idx(1.0, x, y);
                            }
                        }
                        
                    }
                }
                _ => {}
            }
        }        

        fluid.step_simulation(adding_dens, dt);

        
        canvas.set_draw_color(Color::RGBA(0, 0, 0, 255));
        canvas.clear();
        for i in 0..NUM_SQUARES {
            for j in 0..NUM_SQUARES {
                let density = fluid.d.idx( i, j);
                let r = 0;
                let g = (255. * density).clamp(0., 255.) as u8;
                let b = (75. * density).clamp(0., 255.) as u8;
                canvas.set_draw_color(Color::RGBA(r, g, b, 255));
                canvas
                    .fill_rect(Rect::new((i as u32 * GRID_SIZE) as i32, (j as u32 * GRID_SIZE) as i32, GRID_SIZE, GRID_SIZE))
                    .expect("could not fill rect");
                
                
                if show_vec_field {
                    let start_x = (i as u32 * GRID_SIZE + GRID_SIZE/2) as i32;
                    let start_y = (j as u32 * GRID_SIZE + GRID_SIZE/2) as i32;
                    let end_x = start_x + (fluid.u.idx(i, j) * 2.).floor() as i32;
                    let end_y = start_y + (fluid.v.idx( i, j) * 2.).floor() as i32;
                    canvas.set_draw_color(Color::RGBA(122, 122, 122, 255));
                    canvas.draw_line((start_x, start_y), (end_x, end_y)).expect("Could not draw arrow");
                }
                
            }
        }
        canvas.present();
    }
}
