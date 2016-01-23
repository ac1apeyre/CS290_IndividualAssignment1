//Purpose: The engine behind the 3D primitive operations for Mini Assignment 1

//////////////////////////////////////////////
///********    Geometry Functions   *******///
//////////////////////////////////////////////
//This is where you will have to fill in some code

//Purpose: Project vector u onto vector v using the glMatrix library
//Inputs: u (vec3), v (vec3)
//Returns: projv (vec3), the projection of u onto v
function projVector(u, v) {
    var scale = (vec3.dot(u, v)/vec3.dot(v, v));//The scale in front of v is (u dot v) / (v dot v)
    var projv = vec3.create(); //Allocate a vector to hold the output
    vec3.scale(projv, v, scale); //Scale v by the appropriate amount
    return projv; //Return the result
}

//Purpose: To compute the perpendicular projection of a vector u onto a vector v
//Inputs: u (vec3), v (vec3)
//Returns: projperpv (vec3), the projection of u onto v
function projPerpVector(u, v) {
    var projv = projVector(u, v);
    var projperpv = vec3.create();
    vec3.subtract(projperpv, u, projv);
    return projperpv;
}

//Purpose: To compute the angle between the vectors ab and ac
//Inputs: a (vec3), b (vec3), c (vec3)
//Returns: angle (radians - float)
function getAngle(a, b, c) {
    //create vectors from the given points
    //given two points A=(ax,ay,az) and B=(bx,by,bz), the vector *from* A *to* B is (bx-ax,by-ay,bz-az)
    var ab = vec3.create(); //allocate a vector "ab"
    var ac = vec3.create(); //allocate a vector "ac"
    vec3.subtract(ab, b, a); //calculate the vector from point a to point b (ab = b-a)
    vec3.subtract(ac, c, a); //calculate the vector from point a to point c (ac = c-a)
    //find angle between the two vectors
    //cos(theta)= (u dot v) / (|u|*|v|)
    var numerator = vec3.dot(ab, ac); //calculate dot product of vectors ab and ac
    var denominator = vec3.len(ab)*vec3.len(ac); //calculate the product of the vectors' magnitudes
    var quotient = numerator/denominator; //divide dot product by the product of the magnitudes
    var theta = Math.acos(quotient); //calculate the inverse cosine to get angle
    return theta;
}


//Purpose: Given three 3D vertices a, b, and c, compute the area of the triangle
//spanned by them
//Inputs: a (vec3), b (vec3), c (vec3)
//Returns: area (float)
function getTriangleArea(a, b, c) {
    //create vectors from the given points
    //given two points A=(ax,ay,az) and B=(bx,by,bz), the vector *from* A *to* B is (bx-ax,by-ay,bz-az)
    var v1 = vec3.create(); //allocate a vector "v1" (ab)
    var v2 = vec3.create(); //allocate a vector "v2" (ac)
    vec3.subtract(v1, b, a); //calculate the vector from point a to point b (ab = b-a)
    vec3.subtract(v2, c, a); //calculate the vector from point a to point c (ac = c-a)
    
    //area = 0.5 * |v1 x v2|
    var cp = vec3.create(); //allocate a vector "cp" for the cross product of v1 and v2
    vec3.cross(cp, v1, v2); //calculate cross product
    var area = 0.5*vec3.len(cp) //calculate area by halving the magnitude of the cp vector (area of parallelogram formed by v1 and v2)
    return area; 
}

//Purpose: For a plane determined by the points a, b, and c, with the plane
//normal determined by those points in counter-clockwise order using the
//right hand rule, decide whether the point d is above, below, or on the plane
//Inputs: a (vec3), b (vec3), c (vec3)
//Returns: 1 if d is above, -1 if d is below, 0 if d is on
function getAboveOrBelow(a, b, c, d) {
    //create plane vectors using points a, b, c
    //given two points A=(ax,ay,az) and B=(bx,by,bz), the vector *from* A *to* B is (bx-ax,by-ay,bz-az)
    var u1 = vec3.create(); //allocate a vector "u1" (ab)
    var u2 = vec3.create(); //allocate a vector "u2" (ac)
    vec3.subtract(u1, b, a); //calculate the vector from point a to point b (ab = b-a)
    vec3.subtract(u2, c, a); //calculate the vector from point a to point c (ac = c-a)
    
    //calculate plane normal using crossproduct
    var norm = vec3.create(); //allocate a vector for the plane normal
    vec3.cross(norm, u1, u2); //calculate norm using cross product
    
    //create vector from point on plane to point d, arbitrarily A to D
    var ad = vec3.create(); //allocate a vector "ad"
    vec3.subtract(u1, d, a); //calculate the vector from point a to point d (ad = d-a)
    
    //calculate dot product between normal and vector ad
    var dp = vec3.dot(norm, ad); //calculate dot product of normal and ad
    
 	return 9;
    ////determine if d is above, below, or on plane ABC based on sign of the dot product
    //if (dp==0){
    //	return 0; // point is on the plane
    //}
    //else if (dp>0){
    //	return 1; // point is on the same side as the normal vector
    //}
    //else{ // dp<0
    //	return -1; // point is on the opposite side of the normal vector
    //}

}


//Purpose: Given a line segment ab and a line segment cd, compute the intersection
//If they don't intersect, return null
//Inputs: a (vec3), b (vec3), c (vec3), d (vec3)
//Returns: intersection (vec3) or null if no intersection
function getLineSegmentIntersection(a, b, c, d) {
    //TODO: Fill this in for task 4
    return null; //This is a dummy for now.  Fill in with the vec3 
    //representing the intersection if it exists.  Only return null if 
    //no intersection interior to both segments
}

//Purpose: Given three points on a triangle abc, compute the triangle circumcenter
//by intersecting two perpendicular bisectors from two different sides, and
//compute the radius of the circumcircle
//Inputs: a (vec3), b (vec3), c (vec3)
//Returns: On object of the form {circumcenter: vec3, R: float (radius)}
function getTriangleCircumcenter(a, b, c) {
    //TODO: Fill this in for task 5
    return {Circumcenter:vec3.fromValues(0, 0, 0), Radius:0.0};  //This is a dummy
    //for now that shows how to return a JSON object from a function.  Replace
    //vec3.fromValues(0, 0, 0) with the true circumcenter and 0.0 with the 
    //true radius
}

//Purpose: Given four points on a 3D tetrahedron, compute the circumsphere
//by intersecting two perpendicular bisectors from two different triangles,
//and compute the radius of the circumsphere
//Inputs: a (vec3), b (vec3), c (vec3), d (vec3)
//Returns: On object of the form {circumcenter: vec3, R: float (radius)}
function getTetrahedronCircumsphere(a, b, c, d) {
    //EXTRA CREDIT
    return {Circumcenter:vec3.fromValues(0, 0, 0), Radius:0.0};
}

///////////////////////////////////////////////////////////////////
///********           Plotting Utilities                 *******///
///////////////////////////////////////////////////////////////////

//This is code that Chris Tralie has written in to help plot the results
//for help debugging.  Feel free to browse the code to see how plot.ly works
//and ask any questions on the forum

//This is the way I hack the axes to be equal
function getAxesEqual(vs) {
    //Determine the axis ranges
    minval = 0;
    maxval = 0;
    for (var i = 0; i < vs.length; i++) {
        for (var j = 0; j < 3; j++) {
            if (vs[i][j] < minval){ minval = vs[i][j]; }
            if (vs[i][j] > maxval){ maxval = vs[i][j]; }
        }
    }
    return {
    x:{ x: [minval, maxval], y: [0, 0], z: [0, 0],
      mode: 'lines', line: {color: '#000000', width: 1}, type: 'scatter3d', name:'xaxis'
    },
    y:{ x: [0, 0], y: [minval, maxval], z: [0, 0],
      mode: 'lines', line: {color: '#000000', width: 1}, type: 'scatter3d', name:'yaxis'
    },
    z:{ x: [0, 0], y: [0, 0], z: [minval, maxval],
      mode: 'lines', line: {color: '#000000', width: 1}, type: 'scatter3d', name:'zaxis'
    }};
}

function getMousePos(canvas, evt) {
	var rect = canvas.getBoundingClientRect();
	return {
	    X: evt.clientX - rect.left,
	    Y: evt.clientY - rect.top
	};
}

