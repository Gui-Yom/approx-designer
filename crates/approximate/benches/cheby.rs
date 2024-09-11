use approximate::cheby::Chebyshev;
use criterion::{black_box, criterion_group, criterion_main, Criterion};

#[inline(never)]
fn f(x: f32) -> f32 {
    const BETA: f32 = 0.018053968510807;
    const ALPHA: f32 = 1.0 + 5.5 * BETA;
    /// threshold = bt709_oetf(BETA)
    const THRESHOLD: f32 = 0.08124285829863521110029445797874;
    if x >= THRESHOLD {
        ((x + (ALPHA - 1.0)) / ALPHA).powf(1.0 / 0.45)
    } else {
        x / 4.5
    }
}

#[inline(never)]
fn eval(cheby: &Chebyshev<8, f32>, x: f32) -> f32 {
    cheby.eval(x)
}

pub fn criterion_benchmark(c: &mut Criterion) {
    let mut g = c.benchmark_group("f");
    g.bench_function("base", |b| b.iter(|| f(black_box(0.5f32))));

    let cheby = Chebyshev::<16, f32>::fit(0.0, 1.0, f, |x| x.cos());
    g.bench_function("cheby 16", |b| b.iter(|| cheby.eval(black_box(0.5))));
    let cheby = cheby.trunc::<14>();
    g.bench_function("cheby 14", |b| b.iter(|| cheby.eval(black_box(0.5))));
    let cheby = cheby.trunc::<12>();
    g.bench_function("cheby 12", |b| b.iter(|| cheby.eval(black_box(0.5))));
    let cheby = cheby.trunc::<10>();
    g.bench_function("cheby 10", |b| b.iter(|| cheby.eval(black_box(0.5))));
    let cheby = cheby.trunc::<8>();
    g.bench_function("cheby 8", |b| b.iter(|| eval(&cheby, black_box(0.5))));
    let cheby = cheby.trunc::<6>();
    g.bench_function("cheby 6", |b| b.iter(|| cheby.eval(black_box(0.5))));
    let cheby = cheby.trunc::<4>();
    g.bench_function("cheby 4", |b| b.iter(|| cheby.eval(black_box(0.5))));

    g.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
