h = 2;
offset = 0.1;
back = -h-offset;

//fixed
translate([0, 0, back+h/2])
linear_extrude(height = 3*h, center = true)
		import(file = "fixed.svg");

//first piece
//union()  {
//difference(){
//linear_extrude(height = h, center = true)
//		import(file = "piece1.svg");
//
//linear_extrude(height = 2*h, center = true)
//		import(file = "piece1_h1.svg");
//}
//translate([0, 0, back])
//linear_extrude(height = h+offset, center = true)
//		import(file = "piece1_p1.svg");
//
//}


//second piece
//translate([0, 0, back]){
//        difference(){
//    difference(){
//linear_extrude(height = h, center = true)
//		import(file = "piece2.svg");
//
//linear_extrude(height = 2*h, center = true)
//		import(file = "piece2_h1.svg");
//}
//linear_extrude(height = 2*h, center = true)
//		import(file = "piece2_h2.svg");
//
//}
//}


//third piece
//union()  {
//difference(){
//linear_extrude(height = h, center = true)
//		import(file = "piece3.svg");
//
//linear_extrude(height = 2*h, center = true)
//		import(file = "piece3_h1.svg");
//}
//translate([0, 0, back+h+offset])
//linear_extrude(height = 3*h+3*offset, center = true)
//		import(file = "piece3_p1.svg");
//
//}

//fourth piece
//translate([0, 0, -back]){
//        difference(){
//    difference(){
//linear_extrude(height = h, center = true)
//		import(file = "piece4.svg");
//
//linear_extrude(height = 2*h, center = true)
//		import(file = "piece4_h1.svg");
//}
//linear_extrude(height = 2*h, center = true)
//		import(file = "piece4_h2.svg");
//
//}
//}


//figth piece
//union()  {
//linear_extrude(height = h, center = true)
//		import(file = "piece5.svg");
//
//translate([0, 0, h])
//linear_extrude(height = h+offset, center = true)
//		import(file = "piece5_p1.svg");
//
//}