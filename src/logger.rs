use env_logger::{Builder, Env};
use std::io::Write;

pub fn init_logger() {
    let env = Env::default().default_filter_or("info");

    Builder::from_env(env)
        .format(|buf, record| {
            // Timestamp: 24 chars (e.g. [26-07-2026 16:04:59.875])
            let now = chrono::Local::now().format("%d-%m-%Y %H:%M:%S%.3f");

            // Format level string
            let level = match record.level() {
                log::Level::Warn => "WARNING",
                log::Level::Error => "ERROR",
                log::Level::Info => "INFO",
                log::Level::Debug => "DEBUG",
                log::Level::Trace => "TRACE",
            };

            // Module/target name (e.g., 'filter', 'fasta', 'summary')
            let target = record.target().split("::").last().unwrap_or(record.target());

            // Alignment specifiers:
            // {:<7}  -> Left-align level within 7 chars ("WARNING" is 7 chars)
            // {:<10} -> Left-align module target within 10 chars (adjust width as needed)
            writeln!(buf, "[{}] - {:<7} - {:<10} - {}", now, level, target, record.args())
        })
        .init();
}