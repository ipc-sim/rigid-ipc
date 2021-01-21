$fs=0.0001; 

axis_radius = 10;

piston_lenght = 15;
piston_radius = 4;
piston_tickness = 1;
inner_lenght = 6;

spacing = 0.05;

scale([100, 100, 100]){

//rotating part
union(){
    difference(){
        cylinder(1, axis_radius, axis_radius, true);
        cylinder(10, 2, 2, true);
    }
    translate([0, axis_radius*.7, 0.5]){
        cylinder(1, 2, 2);
        translate([0,0,1])
        cylinder(2, 1, 1);
        
        
    translate([0, 0, 2.6])
        cylinder(1, 2, 2);
    }

    translate([0, 0, -4]){
        cylinder(10, 2, 2, true);
   
    }
}


height = 2.3;




//axis
union(){
translate([0, 0, height]){
    translate([0, -axis_radius*.5, 0]){
        cube(size = [1, 2.1*axis_radius, 1], center = true);
        
        translate([0, -2.1*axis_radius/2-1.6/2, 0])rotate([0, 0, 0]){
            difference(){
            cylinder(1.5, 1.5, 1.5, true);
            cylinder(4, 0.7, 0.7, true);
            }
        }
    }
    
    translate([0, axis_radius*.7, 0]){
        difference(){
        cylinder(1.5, 1.8, 1.8, true);
        cylinder(2, 1+spacing, 1+spacing, true);
        }
    }
    

        
}
}


//piston holder
translate([0, 0, 0])
difference(){
rotate([90, 0,0])translate([0, height, axis_radius/2+piston_lenght]){
    difference(){
cylinder(piston_lenght, piston_radius, piston_radius, true);
translate([0, 0, -piston_tickness]){
   cylinder(2*piston_lenght, piston_radius-piston_tickness, piston_radius-piston_tickness, true);
   }
  }
};
translate([0,-12.5,height])rotate([0, 90, 0])translate([0,0,-5])cylinder(10, 1.5, 1.5);
}


//piston
inner_piston_r = piston_radius - piston_tickness - spacing;
difference(){
translate([0, 0, height])rotate([0, 90, 0])translate([0,0,-height])union(){
rotate([90, 0,0])translate([0, height, axis_radius/2+piston_lenght-inner_lenght/2]){
difference(){
    cylinder(inner_lenght, inner_piston_r, inner_piston_r, true);
    translate([0, 0, -piston_tickness]){
        cylinder(inner_lenght, inner_piston_r-piston_tickness, inner_piston_r-piston_tickness, true);
    }
 }


}

translate([0, -axis_radius*.5-2.1*axis_radius/2-1.6/2, height])rotate([0, 90, 0]){
cylinder(2*(inner_piston_r-piston_tickness), 0.7-spacing, 0.7-spacing, true);
}
}

translate([0,-14,height])rotate([0, 90, 0])translate([0,0,-5])cylinder(10, 1.5, 1.5);

}
}