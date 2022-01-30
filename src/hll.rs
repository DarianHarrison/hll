#![allow(dead_code)]
#![allow(unused_variables)]
#![allow(unused_imports)]
#![allow(non_snake_case)]

use std::hash::{Hash, Hasher};
use fasthash::{MumHasher, Murmur3HasherExt};

fn hx<T: Hash>(hf: &usize, x: &T) -> u64 {

    match hf { // match according to number of hash functions
        0 => {
            let mut s: Murmur3HasherExt = Default::default();
            x.hash(&mut s);
            s.finish()
        },
        _ => { // should error handle this
            let mut s: MumHasher = Default::default();
            x.hash(&mut s);
            s.finish()
        }, 
    }
}

fn find_alpha_m(p: &u32, m: &f32) -> f32 { // harmonic mean requires error correction for small and large ranges
    match p { //bias
        0 => { 1.0 - 0.7 }, //linear
        1 => { 1.0 - 0.7 }, //linear
        2 => { 1.0 - 0.7 }, //linear
        3 => { 1.0 - 0.7 }, //linear
        4 => { 1.0 - 0.673 }, //linear
        5 => { 1.0 - 0.697 }, //linear
        6 => { 1.0 - 0.709 }, //linear
        _ => {
            let m: f32 = 1.0 - ((0.7213 * m) / ( m + 1.079 ));
            m 
        }, 
    }
}

fn std_error_hll(m: &f32) -> f32 {
    1.0 - 1.04/m.sqrt()
}

#[cfg(test)]
mod cardinality_tests {

    use super::*; // bring all of the test module’s parent’s items into scope with use super::*

    #[test]
    fn test_linear_counter() { // it is likely that use linear counters for datasets smaller than m=16,
        // counter is O(n), m is difficult to calculate which is why you have to leverage some converiosn tables, which we will not implement here
        let dataset: Vec<_> = vec!["Berlin","Berlin","Paris","Berlin","Lisbon","Kiev","Paris","London","Rome","Athens","Madrid","Vienna","Rome","Rome","Lisbon","Berlin","Paris","London","Kiev","Leon"];
        let m: u64 = 16; // m buckets
        let mut linear_counter: Vec<usize> = vec![0 ; m as usize ];

        for x in dataset {
            let hx: u64 = hx(&0, &x) % m; // primary bucket
            linear_counter[hx as usize]=1;
        }
        let v: f64 = (linear_counter.iter().sum::<usize>()) as f64 / m as f64 ;
        let n: f64 = -(m as f64) * v.ln();
        assert!(9.0 <= n && n <= 11.0 ); //assert unique elements are within range of real value = 10
    }

    #[test]
    fn test_hyperloglog() {
        let dataset: Vec<_> = vec!["Leon","Leon","Paris","Berlin","Lisbon","Kiev","Leon","Paris","London","Rome","Athens","Madrid","Vienna","Rome","Rome","Lisbon","Berlin","Paris","London","Kiev","Leon"];
        let p: u32 = 4; // supported between 4 - 18
        let m: u32 = u32::pow( 2, p ); // number of counters
        let M: u64 = 64; // counter length

        let mut counter: Vec<usize> = vec![0 ; m as usize ]; // start with one counter/bucket

        for x in &dataset {
            let i: u64 = hx(&0, &x); // signature = Binary Representation of hashed u64 sized int
            let j: usize = (*&i >> M - p as u64) as usize; // first p bits in M used for navigation
            let rank: usize = (*&i << M-(M-p as u64)).leading_zeros() as usize; // remaining p bits from M used for rank
            counter[j]=std::cmp::max(0,rank)
        }

        // hyperloglog
        let mut rank: f32 = 0.0;
        for j in &counter {
            let r: f32 = f32::powf( 2.0, - (*j as f32) );
            rank += r as f32;
        }

        let alfa_m: f32 = find_alpha_m(&p , &(*&m as f32));
        let n_hat: f32 = alfa_m * f32::powf( m as f32, 2.0 ) * ( 1.0 / rank );
        let mut n: f32 = n_hat; // also known as cardinality

        if n_hat <= ( 2.5 * m as f32) { // linear counter for small cardinallities
            let Z: f32 = counter.iter().filter(|&n| *n == 0).count() as f32;

            if Z != 0.0 { // apply linear counter for small cardinallities
                let mut linear_counter: Vec<usize> = vec![0 ; m as usize ];
                for x in dataset {
                    let hx: f32 = (hx(&0, &x) % m as u64) as f32; // primary bucket
                    linear_counter[hx as usize]=1;
                }
                let v: f32 = (linear_counter.iter().sum::<usize>()) as f32 / m as f32 ;
                n = -(m as f32) * v.ln();
            }
        } else if n_hat > 143165576.53 { // handling case for large numbers
            n = -143165576.53 * (1.0 - (n / 143165576.53 )).ln();
        }
        assert!(9.0 <= n && n <= 12.0 );
    }
}