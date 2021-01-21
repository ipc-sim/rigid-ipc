
wheel_h = 5;
offset = 0.1;

//curved wheel
//union()  {
//difference(){
//linear_extrude(height = wheel_h, center = true)
//		import(file = "wheel1.svg");
//
//linear_extrude(height = 2*wheel_h, center = true)
//		import(file = "wheel2.svg");
//}
//linear_extrude(height = wheel_h, center = true)
//		import(file = "wheel3.svg");
//
//}



//walls
//translate([0, 0, wheel_h+offset])
//linear_extrude(height = 5*wheel_h, center = true)
//    import(file = "wheel_walls.svg");
    


//piston
//union()  {
//translate([0, 0, -wheel_h - offset])
//linear_extrude(height = wheel_h, center = true)
//		import(file = "wheel_piston1.svg");
//
//translate([0, 0, -offset])
//linear_extrude(height = wheel_h+7*offset, center = true)
//		import(file = "wheel_piston2.svg");
//
//translate([0, 0, wheel_h+offset]){
//
//
//linear_extrude(height = wheel_h, center = true)
//		import(file = "wheel_piston1.svg");
//
//linear_extrude(height = wheel_h, center = true)
//		import(file = "wheel_piston3.svg");
//    
//linear_extrude(height = wheel_h, center = true)
//		import(file = "wheel_piston4.svg");
//    
//
//difference(){
//linear_extrude(height = wheel_h, center = true)
//		import(file = "wheel_piston5.svg");
//    
//linear_extrude(height = 2*wheel_h, center = true)
//		import(file = "wheel_piston6.svg");
//}    
//}
//}
//
//
//
//axis
union(){
translate([0, 0, wheel_h+2*offset]){
    linear_extrude(height = 2*wheel_h, center = true)
		import(file = "axis_1.svg");
    
    translate([0, 0, -wheel_h-2*offset])
    linear_extrude(height = wheel_h, center = true)
		import(file = "wheel_piston5.svg");
    
    translate([0, 0, wheel_h+2*offset])
    linear_extrude(height = wheel_h, center = true)
		import(file = "axis_2.svg");
    
        linear_extrude(height = 2*wheel_h, center = true)
		import(file = "axis_3.svg");
}
}
//
////
////
////
//////end
////union(){
////translate([0, 0, -wheel_h]){
////linear_extrude(height = wheel_h, center = true)
////		import(file = "end_1.svg");
////
////translate([0, 0, wheel_h]){
////linear_extrude(height = wheel_h+2*offset, center = true)
////		import(file = "end_2.svg");
////
////translate([0, 0, wheel_h+offset/2]){
////    difference(){
////linear_extrude(height = wheel_h, center = true)
////		import(file = "end_3.svg");
////linear_extrude(height = 2*wheel_h, center = true)
////		import(file = "end_4.svg");
////}
////}
////}
////}
////}