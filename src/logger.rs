use env_logger::{Builder, Env};
use colored::Colorize;
use std::io::Write;

pub fn init_logger() {
    let env = Env::default().default_filter_or("info");

    Builder::from_env(env)
        .format(|buf, record| {
            let now = chrono::Local::now()
                .format("%Y-%m-%d %H:%M:%S%.3f")
                .to_string()
                .bright_black();

            let (r, g, b) = (106, 171, 115); // Dark Green
            let target = record.target().truecolor(r, g, b);

            let level = match record.level() {
                log::Level::Error => "ERROR  ".truecolor(243, 139, 168),
                log::Level::Warn  => "WARNING".truecolor(255, 165, 0),
                log::Level::Info  => "INFO   ".truecolor(166, 227, 161),
                log::Level::Debug => "DEBUG  ".truecolor(52, 152, 219),
                log::Level::Trace => "TRACE  ".truecolor(155, 89, 182),
            };

            writeln!(buf, "[{}] - {} - {:<25} - {}", now, level, target, record.args())
        })
        .init();
}