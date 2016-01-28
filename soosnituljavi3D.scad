// v mm
Sz = .8;					// presek žice

Nr1 = 2;					// število plasti ovojev
Nz1 = 13;					// število ovojev
Rin1 = 2;					// radij tuljavnika
Zoff1 = 5;				// postavitev tuljavnika

Nr2 = 13;					// število plasti ovojev
Nz2 = 2;					// število ovojev
Rin2 = 3;					// radij tuljavnika
Zoff2 = 0.5;				// postavitev tuljavnika

R = 2;
l = 9;
Zoff = 4.5;

d = 2 * sqrt(Sz / 3.14);
difference() {
	union() {
		tuljava ( Rin1, Zoff1, d, Nz1, Nr1, "green" );
		tuljava ( Rin2, Zoff2, d, Nz2, Nr2, "blue" );
		magnetik ( R, Zoff, l, "red" );
	}
	translate([0, 0, -1]){
	color("white"){
		cube(Nz1*d+Zoff1+2);
	}
	}
}
module magnetik ( R, Zoff, l, barva ) {
	color(barva){
		translate([0, 0, Zoff]){
			cylinder(h = l, r = R, center = true, $fn = 100);
		}
	}
}

module tuljava ( Rin, Zoff, d, Nz, Nr, barva ) {
	Rout = Rin + Nr*d;
	l = Nz*d;
	rad = d/2;
	translate([0, 0, Zoff]){
		for (r = [Rin+rad:d:Rout-rad], z = [-l/2+rad:d:l/2-rad] ){
			zanka ( r, z, rad, barva );
		}
	}
}

module zanka ( R, Z, radij, barva ) {
	color(barva){
		rotate_extrude(angle = 270, convexity = 10, $fn = 100)
		translate([R, Z, 0]){
			circle(r = radij, $fn = 100);
		}
	}
}