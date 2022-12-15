use crate::circuit::base_chip::{BaseChipOps, MUL_COLUMNS, VAR_COLUMNS};
use crate::context::Context;
use crate::pair;
use crate::tests::{random_fr, run_circuit_on_bn256};
use ark_std::{end_timer, start_timer};
use halo2_proofs::pairing::bn256::Fr;

#[test]
fn test_one_line_st() {
    let vars = [(); VAR_COLUMNS].map(|_| random_fr());
    let coeffs = [(); VAR_COLUMNS].map(|_| random_fr());
    let muls_coeffs = [(); MUL_COLUMNS].map(|_| random_fr());
    let next_var = random_fr();
    let next_coeff = random_fr();

    let result: Fr = {
        let mut result = Fr::zero();
        for i in 0..VAR_COLUMNS {
            result = result + vars[i] * coeffs[i]
        }
        for i in 0..MUL_COLUMNS {
            result = result + muls_coeffs[i] * vars[i * 2] * vars[i * 2 + 1]
        }
        result + next_var * next_coeff
    };

    let mut ctx = Context::<Fr>::new();
    {
        let timer = start_timer!(|| "setup");
        for _ in 0..10000 {
            ctx.one_line(
                (0..VAR_COLUMNS)
                    .map(|i| pair!(vars[i], coeffs[i]))
                    .collect(),
                Some(-result),
                (muls_coeffs.to_vec(), Some(next_coeff)),
            );

            ctx.one_line_with_last(vec![], pair!(next_var, Fr::zero()), None, (vec![], None));
        }
        end_timer!(timer);
    }

    run_circuit_on_bn256(ctx, 20);
}

#[test]
fn test_one_line_mt() {
    let vars = [(); VAR_COLUMNS].map(|_| random_fr());
    let coeffs = [(); VAR_COLUMNS].map(|_| random_fr());
    let muls_coeffs = [(); MUL_COLUMNS].map(|_| random_fr());
    let next_var = random_fr();
    let next_coeff = random_fr();

    let result = {
        let mut result = Fr::zero();
        for i in 0..VAR_COLUMNS {
            result = result + vars[i] * coeffs[i]
        }
        for i in 0..MUL_COLUMNS {
            result = result + muls_coeffs[i] * vars[i * 2] * vars[i * 2 + 1]
        }
        result + next_var * next_coeff
    };

    let mut threads = vec![];

    let timer = start_timer!(|| "setup");
    let c = 10000;
    let n = 10;
    let mut ctx = Context::<Fr>::new();
    for i in 0..n {
        let step = c / n;
        let start = i * step * 2;
        let mut ctx = ctx.clone();
        ctx.base_offset = start;
        let t = std::thread::spawn(move || {
            for _ in 0..step {
                ctx.one_line(
                    (0..VAR_COLUMNS)
                        .map(|i| pair!(vars[i], coeffs[i]))
                        .collect(),
                    Some(-result),
                    (muls_coeffs.to_vec(), Some(next_coeff)),
                );

                ctx.one_line_with_last(vec![], pair!(next_var, Fr::zero()), None, (vec![], None));
            }
        });
        threads.push(t);
    }

    for t in threads {
        t.join().unwrap();
    }
    end_timer!(timer);

    ctx.base_offset = c;
    run_circuit_on_bn256(ctx, 20);
}
