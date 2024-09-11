use approximate::cheby::Chebyshev;
use eframe::egui::{CentralPanel, Context, Slider};
use eframe::{Frame, NativeOptions};
use egui_plot::{Legend, Line, PlotPoints};

fn main() {
    eframe::run_native("Approx Designer", NativeOptions::default(), Box::new(|cc| Ok(Box::new(App::new())))).unwrap()
}

struct App {
    cheby: Chebyshev<64, f64>,
    trunc: usize,
    quant: usize,
}

impl App {
    fn new() -> Self {
        let cheby = Chebyshev::<64, _>::fit(0.0, 1.0, Self::eval, |x| x.cos());
        // let cheby = Chebyshev::<16, _>::fit_match(0.0, 1.0, Self::eval, |x| x.cos());
        Self {
            cheby,
            trunc: 8,
            quant: 8,
        }
    }

    #[inline(never)]
    fn eval(x: f64) -> f64 {
        const BETA: f64 = 0.018053968510807;
        const ALPHA: f64 = 1.0 + 5.5 * BETA;
        /// threshold = bt709_oetf(BETA)
        const THRESHOLD: f64 = 0.08124285829863521110029445797874;
        if x >= THRESHOLD {
            ((x + (ALPHA - 1.0)) / ALPHA).powf(1.0 / 0.45)
        } else {
            x / 4.5
        }
    }

    #[inline(never)]
    fn cheby_eval(&self, x: f64) -> f64 {
        self.cheby.eval_trunc(x, self.trunc)
        // self.cheby.const_eval(x)
    }
}

const NUM_POINTS: usize = 8192;

impl eframe::App for App {
    fn update(&mut self, ctx: &Context, frame: &mut Frame) {
        CentralPanel::default().show(ctx, |ui| {
            egui_plot::Plot::new("main").height(ui.available_height() - 40.0).legend(Legend::default()).show(ui, |ui| {
                // ui.set_plot_bounds(PlotBounds::from_min_max([-0.1, -0.1], [1.1, 1.1]));

                let orig_line = Line::new(PlotPoints::from_explicit_callback(Self::eval, 0.0..=1.0, NUM_POINTS)).name("orig");
                let cheby_line = Line::new(PlotPoints::new((0..NUM_POINTS).map(|x| {
                    let x = x as f64 / NUM_POINTS as f64;
                    [x, self.cheby_eval(x)]
                }).collect())).name("cheby");
                let err_line = Line::new(PlotPoints::new((0..NUM_POINTS).map(|x| {
                    let x = x as f64 / NUM_POINTS as f64;
                    let err = if x == 0.0 {
                        Self::eval(x) - self.cheby_eval(x)
                    } else {
                        (Self::eval(x) - self.cheby_eval(x)) / Self::eval(x)
                    };
                    [x, err]
                }).collect())).name("err");
                let quant_orig_line = Line::new(PlotPoints::new((0..NUM_POINTS).map(|x| {
                    let x = x as f64 / NUM_POINTS as f64;
                    [x, quant(Self::eval(x), self.quant)]
                }).collect())).name("quant_orig");
                let quant_cheby_line = Line::new(PlotPoints::new((0..NUM_POINTS).map(|x| {
                    let x = x as f64 / NUM_POINTS as f64;
                    [x, quant(self.cheby_eval(x), self.quant)]
                }).collect())).name("quant_cheby");
                let quant_err_line = Line::new(PlotPoints::new((0..NUM_POINTS).map(|x| {
                    let x = x as f64 / NUM_POINTS as f64;
                    let err = quant(Self::eval(x), self.quant) - quant(self.cheby_eval(x), self.quant);
                    [x, err]
                }).collect())).name("quant_err");
                // let quant_line = Line::new(PlotPoints::new((0..NUM_POINTS).map(|x| {
                //     let x = x as f64 / NUM_POINTS as f64;
                //     let err = quant(x, self.quant);
                //     [x, err]
                // }).collect())).name("quant");
                ui.line(orig_line);
                ui.line(cheby_line);
                ui.line(err_line);
                ui.line(quant_orig_line);
                ui.line(quant_cheby_line);
                ui.line(quant_err_line);
                // ui.line(quant_line);
            });
            ui.add(Slider::new(&mut self.trunc, 1..=64));
            ui.add(Slider::new(&mut self.quant, 8..=16));
        });
    }
}

fn quant(x: f64, bits: usize) -> f64 {
    let scaler = ((1 << bits) - 1) as f64;
    (x * scaler).round() / scaler
}
