/*

Sprocket generator v2

This code is based on the code written by *Talon_1* who based his code on the work of
*Aleksejs*. Big thanks for your contributions. The aim of this code is to be easier
understood by folks that are new to OpenSCAD. The rendered sprocket can be downloaded 
as a .STL file and 3D-printed directly using any slicing program.

*/

//////////////////////
/* CHAIN-PARAMETERS */
//////////////////////

// THESE ARE FOR 25H/04C
roller_d  = 4;
thickness = 1.66;
pitch     = 7;
tolerance = 0.05;

///////////////
/* VARIABLES */
///////////////

teeth    = 8;

// Shaft
bottom_shaft_d = 20;
bottom_shaft_h = 0; // = 0 to remove

top_shaft_d    = 40;
top_shaft_h    = 0; // = 0 to remove

toptop_shaft_d = 36;
toptop_shaft_h = 0; // = 0 to remove

// Bore
hole_d   = 8;

// Holes
number_of_holes = -1;
hole_dia = 4;
hole_ring_dia = 30;

///////////////////////
// RENDERING QUALITY */
///////////////////////

// HIGH : fs=0.25 : fa=3 : fn=0
// LOW  : fs=1    : fa=7 : fn=0
fs = 0.25; // Minimum size of a fragment
fa = 3; // Minimum angle for a fragment
fn = 0; // Number of fragments (overrides fs & fa if non zero)

///////////////
/* MAIN CODE */
///////////////

difference()
{
    // Create a union of four shapes, 3 cylinders and 1 sprocket
    union()
    {
        // Create sprocket using the difined module
        sprocket(teeth, roller_d, pitch, thickness, tolerance);
        
        // Create cylinder on front side of sprocket
        translate([0, 0, thickness])
            cylinder(top_shaft_h, top_shaft_d/2, top_shaft_d/2, $fs=fs, $fa=fa, $fn=fn);
        
        // Create cylinder on back side of sprocket
        rotate([0,180])
            cylinder(bottom_shaft_h, bottom_shaft_d/2, bottom_shaft_d/2, $fs=fs, $fa=fa, $fn=fn);
        
        // Create cylinder on top of the front side cylinder
        translate([0, 0, thickness+top_shaft_h])
            cylinder(toptop_shaft_h, toptop_shaft_d/2, toptop_shaft_d/2, $fs=fs, $fa=fa, $fn=fn);
    }
    
    // Rest of shapes are removal of material
    // Drills out the center hole with 1 mm extra in both directions
    translate([0, 0, -bottom_shaft_h-1])
    {        
        cylinder(bottom_shaft_h+thickness+top_shaft_h+toptop_shaft_h+2, hole_d/2, hole_d/2, $fs=fs, $fa=fa, $fn=fn);
    }

    // Drills 'number_of_holes' many holes in a circle
    angle_between_holes = 360/number_of_holes;
    for(hole_angle = [0:360/number_of_holes:360])
    {
        translate([hole_ring_dia/2*cos(hole_angle), hole_ring_dia/2*sin(hole_angle), -bottom_shaft_h-1])
        {
            cylinder(h = bottom_shaft_h+thickness+top_shaft_h+toptop_shaft_h+2, r = hole_dia/2, $fs=fs, $fa=fa, $fn=fn);
        }
    }
}

/////////////////////
/* SPROCKET MODULE */
/////////////////////

module sprocket(teeth=20, roller=3, pitch=17, thickness=3, tolerance=0.2)
{
	roller_radius = roller/2; //We need radius in our calculations, not diameter
	distance_from_center = pitch/(2*sin(180/teeth));
	angle = (360/teeth);
	
    pitch_radius = sqrt((distance_from_center*distance_from_center) - (pitch*(roller_radius+tolerance))+((roller_radius+tolerance)*(roller_radius+tolerance)));
	    
    difference()
    {
		union()
        {
            // Quality parameters
            $fs = fs; 
            $fa = fa;
            $fn = fn;
            
            // Create inner cylinder with radius = pitch_radius
			cylinder(r=pitch_radius, h=thickness);
            
            // Create outer part of the teeth
			for(tooth=[1:teeth])
            {
				intersection()
                {
					rotate(a=[0, 0, angle*(tooth+0.5)])
                    {
						translate([distance_from_center, 0, 0])
                        {
                            $fs = fs; 
                            $fa = fa;
                            $fn = fn;
							cylinder(r=pitch-roller_radius-tolerance, h=thickness);
						}
					}
					rotate(a=[0,0,angle*(tooth-0.5)])
                    {
						translate([distance_from_center,0,0])
                        {
							cylinder(r=pitch-roller_radius-tolerance,h=thickness);
						}
					}
				}
			}
		}
        
        // Cuts away the inner groove between the teeth
		for(tooth=[1:teeth])
        {
			rotate(a=[0, 0, angle*(tooth+0.5)])
            {
				translate([distance_from_center, 0, -1])
                {
					$fs = fs; 
                    $fa = fa;
                    $fn = fn;
                    cylinder(r=roller_radius+tolerance, h=thickness+2);
				}
			}
		}
	}
}


