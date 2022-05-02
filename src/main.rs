extern crate sdl2;


use sdl2::event::Event;
use sdl2::keyboard::Keycode;
use sdl2::pixels::{Color, PixelFormatEnum};
use sdl2::rect::{Point, Rect};

fn main() {
    let window_width = 1024;
    let grid_size = 32;

    let sdl_context = sdl2::init().unwrap();
    let video_subsystem = sdl_context.video().unwrap();
    let window = video_subsystem
        .window("Fluid Simulation", window_width, window_width)
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
    let creator = canvas.texture_creator();
    let mut texture = creator
        .create_texture_target(PixelFormatEnum::RGBA8888, window_width/grid_size, window_width/grid_size)
        .map_err(|e| e.to_string()).unwrap();

    'mainloop: loop {
        for event in sdl_context.event_pump().unwrap().poll_iter() {
            match event {
                Event::KeyDown {
                    keycode: Some(Keycode::Escape),
                    ..
                } | Event::Quit {..} => break 'mainloop,
                _ => {}
            }
        }
        canvas
            .with_texture_canvas(&mut texture, |texture_canvas| {
                texture_canvas.clear();
                texture_canvas.set_draw_color(Color::RGBA(255, 0, 0, 255));
                texture_canvas
                    .fill_rect(Rect::new(0, 0, grid_size, grid_size))
                    .expect("could not fill rect");
            })
            .map_err(|e| e.to_string()).unwrap();
        canvas.set_draw_color(Color::RGBA(0, 0, 0, 255));
        canvas.clear();
        
        for i in 0..(window_width/grid_size) {
            for j in 0..(window_width/grid_size) {
                canvas.set_draw_color(Color::RGBA((i % 256) as u8, (j % 256) as u8, 0, 255));
                canvas
                    .fill_rect(Rect::new((i * grid_size) as i32, (j * grid_size) as i32, grid_size, grid_size))
                    .expect("could not fill rect");
            }
        }
        
        canvas.present();
    }
}
