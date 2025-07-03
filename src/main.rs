use nalgebra::{Vector3};

use crate::{classes::{molecule::Molecule, primitive_gaussian::PrimitiveGaussian}, energy::nuclear_nuclear_repulsion_energy::nuclear_nuclear_repusion_energy, integrals::{electron_electron_repulsion, electron_nuclear_attraction, kinetic, overlap}};
// use clap::Parser;
mod classes;
mod energy;
mod integrals;
/// Simple program to greet someone
// #[derive(Parser, Debug)]
// #[command(author, version, about, long_about = None)]
// struct Args {

//     /// Name of the person to greet
//     #[arg(short, long, default_value = "Jacob")]
//     name: String,

//     /// Number of times to greet
//     #[arg(short, long, default_value_t = 1)]
//     count: u8,

//     /// DIIS enabled?
//     #[arg(short, long, default_value = "true")]
//     diis: String,
// }

fn main() {
    let sto_3: Vec<(f64, f64)> = vec![
        (0.3425250914E+01, 0.1543289673E+00),
        (0.6239137298E+00, 0.5353281423E+00),
        (0.1688554040E+00, 0.4446345422E+00),
    ];

    let h1_1s: Vec<PrimitiveGaussian> = sto_3.iter()
        .map(|(alpha, coeff)| PrimitiveGaussian::new(*alpha, *coeff, Vector3::new(0., 0., 0.)))
        .collect();

    let h2_1s: Vec<PrimitiveGaussian> = sto_3.iter()
        .map(|(alpha, coeff)| PrimitiveGaussian::new(*alpha, *coeff, Vector3::new(0., 0., 1.)))
        .collect();


    let h1 = classes::atom::Atom::new(1, String::from("H"), h1_1s);
    let h2 = classes::atom::Atom::new(1, String::from("H"), h2_1s);

    let mol = vec![h1, h2];

    let h_mol = Molecule { atoms: mol };

    let e_nn = nuclear_nuclear_repusion_energy(&h_mol);
    let s_mat = overlap::overlap(&h_mol);
    let t_mat = kinetic::kinetic(&h_mol);
    let v_ne = electron_nuclear_attraction::electron_nuclear_attraction(&h_mol);
    let v_ee = electron_electron_repulsion::electron_electron_repulsion(&h_mol);

    // println!("s_mat: \n{}", s_mat);
    // println!("t_mat: \n{}", t_mat);
    // println!("v_ne: \n{}", v_ne);
    println!("v_ee: \n{}", v_ee);
    
    // let args = Args::parse();

    // for _ in 0..args.count {
    //     println!("Hello {}!", args.name);
    // }
    // let diis_flag = args.diis.to_lowercase();

    // if diis_flag == "true" || diis_flag == "t" {
    //     println!("DIIS Enabled!")
    // } else if diis_flag == "false" || diis_flag == "f" {
    //     println!("DIIS Disabled :(")
    // } else {
    //     eprintln!("Unknown command line option for DIIS: {}", args.diis);
    //     std::process::exit(1);
    // }
}
