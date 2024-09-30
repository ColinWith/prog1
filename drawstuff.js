
/* constant globals */

/** number of ray bounces to perform per ray trace */
var NUMBOUNCES = 0; 

var showCustom = false;

/* classes */ 

// Color constructor
class Color {
    constructor(r,g,b,a) {
        try {
            if ((typeof(r) !== "number") || (typeof(g) !== "number") || (typeof(b) !== "number") || (typeof(a) !== "number"))
                throw "color component not a number";
            else if ((r<0) || (g<0) || (b<0) || (a<0)) 
                throw "color component less than 0";
            else if ((r>255) || (g>255) || (b>255) || (a>255)) 
                throw "color component bigger than 255";
            else {
                this.r = r; this.g = g; this.b = b; this.a = a; 
            }
        } // end try
        
        catch (e) {
            console.log(e);
        }
    } // end Color constructor

        // Color change method
    change(r,g,b,a) {
        try {
            if ((typeof(r) !== "number") || (typeof(g) !== "number") || (typeof(b) !== "number") || (typeof(a) !== "number"))
                throw "color component not a number";
            else if ((r<0) || (g<0) || (b<0) || (a<0)) 
                throw "color component less than 0";
            else if ((r>255) || (g>255) || (b>255) || (a>255)) 
                throw "color component bigger than 255";
            else {
                this.r = r; this.g = g; this.b = b; this.a = a; 
            }
        } // end throw
        
        catch (e) {
            console.log(e);
        }
    } // end Color change method

    toString() {
        return (this.r + ", " + this.g + ", " + this.b);
    }
} // end color class

class Vector {
    constructor(x, y, z) {
        this.set(x, y, z);
    }

    set(x, y, z) {
        try {
            if ((typeof (x) !== "number") || (typeof (y) !== "number") || (typeof (z) !== "number")) {
                throw "vector component not a number";
            }
            else {
                this.x = x; this.y = y; this.z = z;
            }
        }

        catch (e) {
            console.log(e);
        }
    }

    copy(v) {
        if (v instanceof Vector) {
            // copy values from Vector3D
            this.x = v.x; this.y = v.y; this.z = v.z;
        }
        else {
            // copy values from array of values
            this.x = v[0]; this.y = v[1]; this.z = v[2];
        }
    }

    // vector math

    static add(v1, v2) {
        try {
            if (v1 instanceof Vector && v2 instanceof Vector) {
                return (new Vector(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z));
            }
            else {
                throw "Vector.add component not a Vector";
            }
        }
        catch (e) {
            console.log(e);
            return (new Vector(NaN, NaN, NaN));
        }
    }
    static sub(v1, v2) {
        try {
            if (v1 instanceof Vector && v2 instanceof Vector) {
                return (new Vector(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z));
            }
            else {
                throw "Vector.sub component not a Vector";
            }
        }
        catch (e) {
            console.log(e);
            return (new Vector(NaN, NaN, NaN));
        }
    }
    static multConst(v, c) {
        try {
            if (v instanceof Vector) {
                return (new Vector(v.x * c, v.y * c, v.z * c));
            }
            else {
                throw "Vector.multConst component not a Vector";
            }
        }
        catch (e) {
            console.log(e);
            return (new Vector(NaN, NaN, NaN));
        }
    }
    static divConst(v, c) {
        try {
            if (v instanceof Vector) {
                return (new Vector(v.x / c, v.y / c, v.z / c));
            }
            else {
                throw "Vector.divConst component not a Vector";
            }
        }
        catch (e) {
            console.log(e);
            return (new Vector(NaN, NaN, NaN));
        }
    }
    // finds the dot product between this vector and the given one
    static dot(v1, v2) {
        try {
            if (v1 instanceof Vector && v2 instanceof Vector) {
                return (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
            }
            else {
                throw "Vector.dot component not a Vector";
            }
        }
        catch (e) {
            console.log(e);
            return null;
        }
    }
    // finds cross product between the given vectors
    static cross(v1, v2) {
        try {
            if (v1 instanceof Vector && v2 instanceof Vector) {
                return (new Vector((v1.y * v2.z) - (v1.z * v2.y),
                    (v1.z * v2.x) - (v1.x * v2.z),
                    (v1.x * v2.y) - (v1.y * v2.x)));
            }
            else {
                throw "Vector.cross component not a Vector";
            }
        }
        catch (e) {
            console.log(e);
            return (new Vector(NaN, NaN, NaN));
        }
    }
    // returns the magnitude of the given vector
    static magnitude(v) {
        try {
            if (v instanceof Vector) {
                return Math.sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
            }
            else {
                throw "Vector.magnitude component not a Vector";
            }
        }
        catch (e) {
            console.log(e);
            return null;
        }
    }
    // returns a normalized version of the given vector
    static normalize(v) {
        try {
            if (v instanceof Vector) {
                return Vector.divConst(v, Vector.magnitude(v));
            }
            else {
                throw "Vector.normalize component not a Vector";
            }
        }
        catch (e) {
            console.log(e);
            return (new Vector(NaN, NaN, NaN));
        }
    }

    // returns true if the given vectors are equal
    static equal(v1, v2) {
        try {
            if (v1 instanceof Vector && v2 instanceof Vector) {
                if ( v1.x == v2.x && v1.y == v2.y && v1.z == v2.z ) {
                    return true;
                }
                else {
                    return false;
                }
            }
            else {
                throw "Vector.equal component not a Vector";
            }
        }
        catch (e) {
            console.log(e);
            return null;
        }
    }

    // returns this vector as a string
    toString() {
        return (this.x + ", " + this.y + ", " + this.z);
    }
}

/* utility functions */

// draw a pixel at x,y using color
function drawPixel(imagedata,x,y,color) {
    try {
        if ((typeof(x) !== "number") || (typeof(y) !== "number"))
            throw "drawpixel location not a number";
        else if ((x<0) || (y<0) || (x>=imagedata.width) || (y>=imagedata.height))
            throw "drawpixel location outside of image";
        else if (color instanceof Color) {
            var pixelindex = (y*imagedata.width + x) * 4;
            imagedata.data[pixelindex] = color.r;
            imagedata.data[pixelindex+1] = color.g;
            imagedata.data[pixelindex+2] = color.b;
            imagedata.data[pixelindex+3] = color.a;
        } else 
            throw "drawpixel color is not a Color";
    } // end try
    
    catch(e) {
        console.log(e);
    }
} // end drawPixel
    
// draw random pixels
function drawRandPixels(context) {
    var c = new Color(0,0,0,0); // the color at the pixel: black
    var w = context.canvas.width;
    var h = context.canvas.height;
    var imagedata = context.createImageData(w,h);
    const PIXEL_DENSITY = 0.01;
    var numPixels = (w*h)*PIXEL_DENSITY; 
    
    // Loop over 1% of the pixels in the image
    for (var x=0; x<numPixels; x++) {
        c.change(Math.random()*255,Math.random()*255,
            Math.random()*255,255); // rand color
        drawPixel(imagedata,
            Math.floor(Math.random()*w),
            Math.floor(Math.random()*h),
                c);
    } // end for x
    context.putImageData(imagedata, 0, 0);
} // end draw random pixels

// get the input ellipsoids from the standard class URL
function getInputEllipsoids() {
    const INPUT_ELLIPSOIDS_URL = 
        "https://ncsucgclass.github.io/prog1/ellipsoids.json";
        
    // load the ellipsoids file
    var httpReq = new XMLHttpRequest(); // a new http request
    httpReq.open("GET",INPUT_ELLIPSOIDS_URL,false); // init the request
    httpReq.send(null); // send the request
    var startTime = Date.now();
    while ((httpReq.status !== 200) && (httpReq.readyState !== XMLHttpRequest.DONE)) {
        if ((Date.now()-startTime) > 3000)
            break;
    } // until its loaded or we time out after three seconds
    if ((httpReq.status !== 200) || (httpReq.readyState !== XMLHttpRequest.DONE)) {
        console.log*("Unable to open input ellipses file!");
        return String.null;
    } else
        return JSON.parse(httpReq.response); 
} // end get input ellipsoids

//get the input triangles from the standard class URL
function getInputTriangles() {
    const INPUT_TRIANGLES_URL = 
        "https://ncsucgclass.github.io/prog1/triangles2.json";
        
    // load the triangles file
    var httpReq = new XMLHttpRequest(); // a new http request
    httpReq.open("GET",INPUT_TRIANGLES_URL,false); // init the request
    httpReq.send(null); // send the request
    var startTime = Date.now();
    while ((httpReq.status !== 200) && (httpReq.readyState !== XMLHttpRequest.DONE)) {
        if ((Date.now()-startTime) > 3000)
            break;
    } // until its loaded or we time out after three seconds
    if ((httpReq.status !== 200) || (httpReq.readyState !== XMLHttpRequest.DONE)) {
        console.log*("Unable to open input triangles file!");
        return String.null;
    } else
        return JSON.parse(httpReq.response); 
} // end get input triangles

//get the input boxex from the standard class URL
function getInputBoxes() {
    const INPUT_BOXES_URL = 
        "https://ncsucgclass.github.io/prog1/boxes.json";
        
    // load the boxes file
    var httpReq = new XMLHttpRequest(); // a new http request
    httpReq.open("GET",INPUT_BOXES_URL,false); // init the request
    httpReq.send(null); // send the request
    var startTime = Date.now();
    while ((httpReq.status !== 200) && (httpReq.readyState !== XMLHttpRequest.DONE)) {
        if ((Date.now()-startTime) > 3000)
            break;
    } // until its loaded or we time out after three seconds
    if ((httpReq.status !== 200) || (httpReq.readyState !== XMLHttpRequest.DONE)) {
        console.log*("Unable to open input boxes file!");
        return String.null;
    } else
        return JSON.parse(httpReq.response); 
} // end get input boxes

// put random points in the ellipsoids from the class github
function drawRandPixelsInInputEllipsoids(context) {
    var inputEllipsoids = getInputEllipsoids();
    var w = context.canvas.width;
    var h = context.canvas.height;
    var imagedata = context.createImageData(w,h);
    const PIXEL_DENSITY = 0.1;
    var numCanvasPixels = (w*h)*PIXEL_DENSITY; 
    
    if (inputEllipsoids != String.null) { 
        var x = 0; var y = 0; // pixel coord init
        var cx = 0; var cy = 0; // init center x and y coord
        var ellipsoidXRadius = 0; // init ellipsoid x radius
        var ellipsoidYRadius = 0; // init ellipsoid y radius
        var numEllipsoidPixels = 0; // init num pixels in ellipsoid
        var c = new Color(0,0,0,0); // init the ellipsoid color
        var n = inputEllipsoids.length; // the number of input ellipsoids
        //console.log("number of ellipses: " + n);

        // Loop over the ellipsoids, draw rand pixels in each
        for (var e=0; e<n; e++) {
            cx = w*inputEllipsoids[e].x; // ellipsoid center x
            cy = h*inputEllipsoids[e].y; // ellipsoid center y
            ellipsoidXRadius = Math.round(w*inputEllipsoids[e].a); // x radius
            ellipsoidYRadius = Math.round(h*inputEllipsoids[e].b); // y radius
            numEllipsoidPixels = ellipsoidXRadius*ellipsoidYRadius*Math.PI; // projected ellipsoid area
            numEllipsoidPixels *= PIXEL_DENSITY; // percentage of ellipsoid area to render to pixels
            numEllipsoidPixels = Math.round(numEllipsoidPixels);
            //console.log("ellipsoid x radius: "+ellipsoidXRadius);
            //console.log("ellipsoid y radius: "+ellipsoidYRadius);
            //console.log("num ellipsoid pixels: "+numEllipsoidPixels);
            c.change(
                inputEllipsoids[e].diffuse[0]*255,
                inputEllipsoids[e].diffuse[1]*255,
                inputEllipsoids[e].diffuse[2]*255,
                255); // ellipsoid diffuse color
            for (var p=0; p<numEllipsoidPixels; p++) {
                do {
                    x = Math.random()*2 - 1; // in unit square 
                    y = Math.random()*2 - 1; // in unit square
                } while (Math.sqrt(x*x + y*y) > 1) // a circle is also an ellipse
                drawPixel(imagedata,
                    cx+Math.round(x*ellipsoidXRadius),
                    cy+Math.round(y*ellipsoidYRadius),c);
                //console.log("color: ("+c.r+","+c.g+","+c.b+")");
                //console.log("x: "+Math.round(w*inputEllipsoids[e].x));
                //console.log("y: "+Math.round(h*inputEllipsoids[e].y));
            } // end for pixels in ellipsoid
        } // end for ellipsoids
        context.putImageData(imagedata, 0, 0);
    } // end if ellipsoids found
} // end draw rand pixels in input ellipsoids

// draw 2d projections read from the JSON file at class github
function drawInputEllipsoidsUsingArcs(context) {
    var inputEllipsoids = getInputEllipsoids();
    
    
    if (inputEllipsoids != String.null) { 
        var c = new Color(0,0,0,0); // the color at the pixel: black
        var w = context.canvas.width;
        var h = context.canvas.height;
        var n = inputEllipsoids.length; 
        //console.log("number of ellipsoids: " + n);

        // Loop over the ellipsoids, draw each in 2d
        for (var e=0; e<n; e++) {
            context.fillStyle = 
                "rgb(" + Math.floor(inputEllipsoids[e].diffuse[0]*255)
                +","+ Math.floor(inputEllipsoids[e].diffuse[1]*255)
                +","+ Math.floor(inputEllipsoids[e].diffuse[2]*255) +")"; // diffuse color
            context.save(); // remember previous (non-) scale
            context.scale(1, inputEllipsoids[e].b/inputEllipsoids[e].a); // scale by ellipsoid ratio 
            context.beginPath();
            context.arc(
                Math.round(w*inputEllipsoids[e].x),
                Math.round(h*inputEllipsoids[e].y),
                Math.round(w*inputEllipsoids[e].a),
                0,2*Math.PI);
            context.restore(); // undo scale before fill so stroke width unscaled
            context.fill();
            //console.log(context.fillStyle);
            //console.log("x: "+Math.round(w*inputEllipsoids[e].x));
            //console.log("y: "+Math.round(h*inputEllipsoids[e].y));
            //console.log("a: "+Math.round(w*inputEllipsoids[e].a));
            //console.log("b: "+Math.round(h*inputEllipsoids[e].b));
        } // end for ellipsoids
    } // end if ellipsoids found
} // end draw input ellipsoids

//put random points in the triangles from the class github
function drawRandPixelsInInputTriangles(context) {
    var inputTriangles = getInputTriangles();
    var w = context.canvas.width;
    var h = context.canvas.height;
    var imagedata = context.createImageData(w,h);
    const PIXEL_DENSITY = 0.1;
    var numCanvasPixels = (w*h)*PIXEL_DENSITY; 
    
    if (inputTriangles != String.null) { 
        var x = 0; var y = 0; // pixel coord init
        var cx = 0; var cy = 0; // init center x and y coord
        var numTrianglePixels = 0; // init num pixels in triangle
        var c = new Color(0,0,0,0); // init the triangle color
        var n = inputTriangles.length; // the number of input files
        //console.log("number of files: " + n);

        // Loop over the triangles, draw rand pixels in each
        for (var f=0; f<n; f++) {
        	var tn = inputTriangles[f].triangles.length;
        	//console.log("number of triangles in this files: " + tn);
        	
        	// Loop over the triangles, draw each in 2d
        	for(var t=0; t<tn; t++){
        		var vertex1 = inputTriangles[f].triangles[t][0];
        		var vertex2 = inputTriangles[f].triangles[t][1];
        		var vertex3 = inputTriangles[f].triangles[t][2];

        		var vertexPos1 = inputTriangles[f].vertices[vertex1];
        		var vertexPos2 = inputTriangles[f].vertices[vertex2];
        		var vertexPos3 = inputTriangles[f].vertices[vertex3];
        		//console.log("vertexPos1 " + vertexPos1);
        		//console.log("vertexPos2 " + vertexPos2);
        		//console.log("vertexPos3 " + vertexPos3);
        		
        		// triangle position on canvas
        		
        		var v1 = [w*vertexPos1[0], h*vertexPos1[1]];
        		var v2 = [w*vertexPos2[0], h*vertexPos2[1]];
        		var v3 = [w*vertexPos3[0], h*vertexPos3[1]];
        		
        		// calculate triangle area on canvas (shoelace formula)
        		var triangleArea = 0.5*Math.abs(v1[0]*v2[1]+v2[0]*v3[1]+v3[0]*v1[1]-v2[0]*v1[1]-v3[0]*v2[1]-v1[0]*v3[1]);
        		var numTrianglePixels = triangleArea; // init num pixels in triangle
            	//console.log("triangle area " + triangleArea);
            	numTrianglePixels *= PIXEL_DENSITY; // percentage of triangle area to render to pixels
            	numTrianglePixels = Math.round(numTrianglePixels);
            	// console.log("numTrianglePixels " + numTrianglePixels);
            	c.change(
            		inputTriangles[f].material.diffuse[0]*255,
                	inputTriangles[f].material.diffuse[1]*255,
                	inputTriangles[f].material.diffuse[2]*255,
                	255); // triangle diffuse color
            	for (var p=0; p<numTrianglePixels; p++) {
                    var point; // on canvas plane
            		var triangleTest = 0;
            		while (triangleTest == 0 ){ //if the pixel outside the triangle
                  
            			point = [Math.floor(Math.random()*w), Math.floor(Math.random()*h)];
                    	// plane checking
            			
                    	var t1 = ((point[0]-v2[0]) * (v1[1] - v2[1]) - (v1[0] - v2[0]) * (point[1] - v2[1])) < 0.0;
                    	var t2 = ((point[0]-v3[0]) * (v2[1] - v3[1]) - (v2[0] - v3[0]) * (point[1] - v3[1])) < 0.0;
                    	var t3 = ((point[0]-v1[0]) * (v3[1] - v1[1]) - (v3[0] - v1[0]) * (point[1] - v1[1])) < 0.0;
                    	
                    	if((t1==t2)&&(t2==t3)) // draw the pixel if inside the triangle
                    		triangleTest = 1;
            		}
            		drawPixel(imagedata,point[0],point[1],c);
                	//console.log("color: ("+c.r+","+c.g+","+c.b+")");
                	//console.log("x: "+ x);
                	//console.log("y: "+ y);
            	} // end for pixels in triangle
        	} // end for triangles
    	} // end for files
        context.putImageData(imagedata, 0, 0);
    } // end if triangle file found
} // end draw rand pixels in input triangles

/**
 * Calculates the ambient vector from the given values
 * @param Ka reflectivity coefficient
 * @param La color
 * @returns ambient vector
 */
function ambient(Ka, La) {
    return Ka * La;
}
/** 
 * Calculates the Diffuse vector from the given values
 * @param Kd reflectivity coefficient
 * @param Ld color
 * @param N normal vector
 * @param L light vector
 * @returns diffuse vector
 */
function diffuse(Kd, Ld, N, L) {
    try {
        if (N instanceof Vector && L instanceof Vector) {
            return Kd * Ld * Math.max(Vector.dot(Vector.normalize(N), Vector.normalize(L)), 0);
        }
        else {
            throw "Diffuse component not a Vector";
        }
    }
    catch (e) {
        console.log(e);
        return NaN;
    }
    
}
/**
 * Calculates the Spectular vector from the given values
 * @param Ks reflectivity coefficient
 * @param Ls color
 * @param N normal vector
 * @param L light vector
 * @param V view point vector
 * @param a reflectivity exponent
 * @returns specular vector
 */
function spectular(Ks, Ls, N, V, L, a) {
    try {
        if (N instanceof Vector && V instanceof Vector && L instanceof Vector) {
            var H = halfVector(V, L);
            return Ks * Ls * Math.pow(Math.max(Vector.dot(Vector.normalize(N), Vector.normalize(H)), 0), a);
        }
        else {
            throw "spectular component not a Vector";
        }
    } catch (e) {
        console.log(e);
        return NaN;
    }
}
/** 
 * Calculates the Reflection vector from the given values
 * @param L light direction vector
 * @param N normal vector
 * @returns reflection vector
 */
function reflection(L, N) {
    try {
        if (L instanceof Vector && N instanceof Vector) {
            return Vector.normalize(Vector.sub(Vector.multConst(Vector.normalize(N), 2.0 * Vector.dot(Vector.normalize(L), Vector.normalize(N))), Vector.normalize(L)));
        }
        else {
            throw "reflection component not a Vector";
        }
    }
    catch (e) {
        console.log(e);
        return new Vector(NaN, NaN, NaN);
    }
}
/** 
 * Calculates the half vector from the given values
 * @param V vision vector
 * @param L light direction vector
 * @returns halfVector
 */
function halfVector(V, L) {
    try {
        if (V instanceof Vector && L instanceof Vector) {
            return Vector.normalize(Vector.add(Vector.normalize(V), Vector.normalize(L)));
        }
        else {
            throw "halfVector component not a Vector";
        }
    }
    catch (e) {
        console.log(e);
        return new Vector(NaN, NaN, NaN);
    }
}

/**
 * Calculates the color of a vertex from the given 
 * parameters, using the blinn-phong model.
 * @param {any} Ka ambient reflectivity coefficient
 * @param {any} La ambient light color
 * @param {any} S shadow ray occluded
 * @param {any} Kd diffuse reflectivity coefficient
 * @param {any} Ld diffuse light color
 * @param {any} N normal
 * @param {any} L light vector
 * @param {any} Ks specular reflectivity coefficient
 * @param {any} Ls specular light color
 * @param {any} V vision vector
 * @param {any} a reflectivity exponent
 * @return vertex color
 */
function vertexColor(Ka, La, S, Kd, Ld, N, L, Ks, Ls, V, a) {
    try {
        if (N instanceof Vector && L instanceof Vector && V instanceof Vector) {
            var ambi = ambient(Ka, La);
            //console.log("ambient: " + ambi);
            var diff = diffuse(Kd, Ld, N, L);
            //console.log("diffuse: " + diff);
            var spec = spectular(Ks, Ls, N, V, L, a);
            //console.log("spectular: " + spec);

            //var col = ambi + diff + spec;
            var col = ambi + (S * (diff + spec));

            //console.log("color: " + col);
            return col;
        }
        else {
            throw "vertexColor component not a Vector";
        }
    }
    catch (e) {
        console.log(e);
        return null;
    }
}

/**
 * Checks if the given point is within the tiangle made up from the given vertices
 * @param {any} point point to check
 * @param {any} v1 first vertex of the triangle
 * @param {any} v2 second vertex of the triangle
 * @param {any} v3 third vertex of the triangle
 * @returns true if given point is within the given triangles vertices
 */
function isInTriangle(point, v1, v2, v3) {
    try {
        if (point instanceof Vector && v1 instanceof Vector && v2 instanceof Vector && v3 instanceof Vector) {
            // plane checking
            var t1 = ((point.x - v2.x) * (v1.y - v2.y) - (v1.x - v2.x) * (point.y - v2.y)) < 0.0;
            var t2 = ((point.x - v3.x) * (v2.y - v3.y) - (v2.x - v3.x) * (point.y - v3.y)) < 0.0;
            var t3 = ((point.x - v1.x) * (v3.y - v1.y) - (v3.x - v1.x) * (point.y - v1.y)) < 0.0;


            if ((t1 == t2) && (t2 == t3)) { // return true to draw the pixel if inside the triangle
                return true;
            }
            return false;
        }
        else {
            throw "isInTriangle component not a Vector";
        }

    }
    catch (e) {
        console.log(e);
        return false; 
    }

}

/**
 * Finds the normal vector from the three given points.
 * @param {any} A first vertex point
 * @param {any} B second vertex point
 * @param {any} C third vertex point
 * @returns normal vector
 */
function triangleNormal(A, B, C) {
    try {
        if (A instanceof Vector && B instanceof Vector && C instanceof Vector) {
            return Vector.cross(Vector.sub(B, A), Vector.sub(C, A));
        }
        else {
            throw "triangleNormal component not a Vector";
        }
    }
    catch(e) {
        console.log(e);
        return new Vector(NaN, NaN, NaN);
    }

}

/**
 * Takes eye position (E), point (P), and three triangle vectors (v1, v2, v3)
 * and check for  an intersection
 * returns point of intersection if true
 * @param {any} start start positon vector
 * @param {any} dir direction vector
 * @param {any} v1 triangle vertex 1
 * @param {any} v2 triangle vertex 2
 * @param {any} v3 triangle vertex 3
 * @returns point of intersection, if any
 */
function rayIntersection(start, dir, v1, v2, v3) {
    try {
        if (start instanceof Vector && dir instanceof Vector && v1 instanceof Vector && v2 instanceof Vector && v3 instanceof Vector) {
            //console.log("E: " + E + ", P: " + P + ", v1: " + v1 + ", v2: " + v2 + ", v3: " + v3);
            //var D = Vector.sub(P, E); // ray direction
            //console.log("D: " + D);
            var N = triangleNormal(v1, v2, v3); // triangle normal
            //console.log("N: " + N);
            var pd = Vector.dot(N, v1); // triangle plane coefficient
            //console.log("pd: " + pd);

            var ND = Vector.dot(N, dir);
            //console.log("ND: " + ND);
            if (ND == 0) {
                return (new Vector(NaN, NaN, NaN)); // no intersection
            }
            else {
                // does intersect
                var t = (pd - Vector.dot(N, start)) / ND; // ray distance to intersection
                //console.log("t: " + t);

                // calculate coordinates of intersection
                var I = Vector.add(start, Vector.multConst(dir, t));
                //console.log("[INTERSECTION AT] I: " + I[0] + ", " + I[1] + ", " + I[2]);
                //return { isectT: t, vector: I };
                return I;

            }
        }
        else {
            throw "rayIntersection component not a Vector";
        }
    }
    catch (e) {
        console.log(e);
        return (new Vector(NaN, NaN, NaN));
    }
}

/**
 * Checks for any intersections from the start point, traveling in the given direction.
 * returns closest surface intersected with, formatted as:
 * { isectPos: , isectVec: , normal: , material: , context:  }
 * 
 * @param {any} inputTriangles
 * @param {any} start start positon vector
 * @param {any} dir direction vector
 * @returns closest surface in intersection
 */
function surfaceRayCast(inputTriangles, start, dir, context) {
    try {
        if (inputTriangles == String.null || !(start instanceof Vector) || !(dir instanceof Vector)) {
            throw "surfaceRayCast component is invalid"
        }
    }
    catch(e) {
        console.log(e);
        return null;
    }

    // TEMP VARIABLE
    //var testAgainstGroupings = true;


    var w = context.canvas.width;
    var h = context.canvas.height;

    var vertex1;
    var vertex2;
    var vertex3;

    var vertexPos1 = new Vector(NaN, NaN, NaN);
    var vertexPos2 = new Vector(NaN, NaN, NaN);
    var vertexPos3 = new Vector(NaN, NaN, NaN);
    // triangle position on canvas

    var v1;
    var v2;
    var v3;

    var newIntersection = null;
    var prevIntersection = null;
    var surface = null;

    var n = inputTriangles.length;

    // Loop over the input files
    for (var f = 0; f < n; f++) {

        var tn = inputTriangles[f].triangles.length;
        //console.log("number of triangles in this files: " + tn);

        // Check each triangle for intersections
        for (var t = 0; t < tn; t++) {
            vertex1 = inputTriangles[f].triangles[t][0];
            vertex2 = inputTriangles[f].triangles[t][1];
            vertex3 = inputTriangles[f].triangles[t][2];
            //console.log("vertexs: (" + vertex1 + "), (" + vertex2 + "), (" + vertex3 + ")");

            vertexPos1.copy(inputTriangles[f].vertices[vertex1]);
            vertexPos2.copy(inputTriangles[f].vertices[vertex2]);
            vertexPos3.copy(inputTriangles[f].vertices[vertex3]);

            //console.log("vertexPos: (" + vertexPos1.toString() + "), (" + vertexPos2.toString() + "), (" + vertexPos3.toString() + ")");

            // triangle position on canvas
            v1 = new Vector(w * vertexPos1.x, h * vertexPos1.y, 0);
            v2 = new Vector(w * vertexPos2.x, h * vertexPos2.y, 0);
            v3 = new Vector(w * vertexPos3.x, h * vertexPos3.y, 0);
            //console.log("v's: (" + v1.toString() + "), (" + v2.toString() + "), (" + v3.toString() + ")");


            //var newI = rayIntersection([eyePos[0], eyePos[1], eyePos[2]], point, vertexPos1, vertexPos2, vertexPos3);
            newIntersection = rayIntersection(start, dir, vertexPos1, vertexPos2, vertexPos3);

            // check if valid
            if (newIntersection != null && !(newIntersection.z < 0)
                && isInTriangle(new Vector(newIntersection.x * w, newIntersection.y * h, 0), v1, v2, v3)
                && !(Vector.equal(newIntersection, start))) {
                // check if closest
                if (prevIntersection == null || newIntersection.z < prevIntersection.z) {
                    surface = { isectPos: newIntersection, isectVec: dir, normal: triangleNormal(vertexPos1, vertexPos2, vertexPos3), material: inputTriangles[f].material , context: context};
                }
                prevIntersection = newIntersection;

            }

        } // end checking shapes
    }
    return surface;
}

function shadowCheck(L, N) {
    var difference = Vector.magnitude(Vector.sub(L, N));
    if ((L.x + L.y + L.z) < 0) {
        difference *= -1;
    }
    if ((N.x + N.y + N.z) < 0) {
        difference *= -1;
    }


    if ( difference > 0 ) {
        
        return 1;
    }
    else {
        //console.log();
        return 0;
    }
}

/**
 * Calculate the color of the given surface
 * @param {any} surface { isectPos: , isectVec: , normal: , material:  , context}
 * @param {any} whichBounce numbered bounce
 * @param {any} viewPoint 
 * @param {any} light 
 * @returns shaded color
 */
function shade(surface, whichBounce, viewPoint, light) {
    try {
        if (surface == null || whichBounce == null || !(viewPoint instanceof Vector) || !(light.pos instanceof Vector)) {
            throw "shade component is invalid"
        }
    }
    catch (e) {
        console.log(e);
        return null;
    }

    var normalVec = surface.normal;
    var lightVec = Vector.sub(light.pos, surface.isectPos);
    var visionVec = Vector.sub(viewPoint, surface.isectPos);
    var shadedColor = [0, 0, 0]; // the shaded color at the pixel
    var surfaceR = null;

    // Calculate shading

    shadowOcc = shadowCheck(lightVec, normalVec);

    // get the vertex color values for this vertex
    for (var i = 0; i < 3; i++) {
        shadedColor[i] = vertexColor(surface.material.ambient[i], light.color[i], shadowOcc,
            surface.material.diffuse[i], light.color[i], normalVec, lightVec,
            surface.material.specular[i], light.color[i], visionVec, surface.material.n);
    }
    shadedColor[3] = 1; // sets alpha value to 100%

    // recursive reflection bounce check
    if (whichBounce < NUMBOUNCES) {
        
        surfaceR = surfaceRayCast(getInputTriangles(), surface.isectPos, reflection(lightVec, surface.normal), surface.context);

        if (surfaceR != null) {
            var shaded = shade(surfaceR, whichBounce + 1, surface.isectPos, light);

            // apply reflection to each of the colors values
            for (var i = 0; i < 3; i++) {
                shadedColor[i] += surface.material.specular[i] * shaded[i];
            }
        }
        
    }

    return shadedColor;

}

// draws the input triangles using pixels
function drawInputTraingles(context) {
    var inputTriangles = getInputTriangles();
    var inputTriangles_custom = [
        // pink
        {
            "material": { "ambient": [0.945, 0.706, 0.455], "diffuse": [0.3, 0.3, 0.3], "specular": [0.15, 0.15, 0.15], "n": 6 },
            "vertices": [
                [-2.25, -2.25, 0.5], [-2.25, 4.25, 0.5], [5.25, 3.0, 0.5]
               

            ],
            "triangles": [[0, 1, 2]]
        },
        // dark pink
        {
            "material": { "ambient": [0.851, 0.514, 0.349], "diffuse": [0.3, 0.3, 0.3], "specular": [0.15, 0.15, 0.15], "n": 10 },
            "vertices": [
                [-1, 0.12, 0.5], [1.5, 0.22, 0.5], [1.0, -1.5, 0.0],

            ],
            "triangles": [[0, 1, 2]]
        },
        // darker pink
        {
            "material": { "ambient": [0.780, 0.353, 0.263], "diffuse": [0.3, 0.3, 0.3], "specular": [0.1, 0.1, 0.1], "n": 6 },
            "vertices": [
                [0.53, 0.29, 0.29], [0.85, 0.29, 0.29], [0.68, 0.25, 0.29],
                [0.52, 0.13, 0.29], [0.53, 0.29, 0.29], [0.58, 0.30, 0.29],
                [0.65, 0.25, 0.29], [0.70, 0.26, 0.29], [0.68, 0.096, 0.29],
                [0.70, 0.29, 0.29], [0.74, 0.30, 0.29], [0.75, 0.13, 0.29],
                [0.80, 0.29, 0.29], [0.85, 0.29, 0.29], [0.88, 0.12, 0.29],
                [0.73, 0.29, 0.29], [0.73, 0.58, 0.29], [0.90, 0.58, 0.29], [0.85, 0.29, 0.29]

            ],
            "triangles": [[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11], [12, 13, 14], [15, 16, 17], [15, 17, 18]]
        },
        // white 1
        {
            "material": { "ambient": [0.941, 0.808, 0.412], "diffuse": [0.3, 0.3, 0.3], "specular": [0.25, 0.25, 0.25], "n": 5 },
            "vertices": [
                [0.43, 0.43, 0.26], [0.53, 0.61, 0.25], [0.61, 0.54, 0.25],
                [0.43, 0.43, 0.26], [0.62, 0.56, 0.28], [0.65, 0.48, 0.28],
                [0.58, 0.54, 0.28], [0.77, 0.63, 0.28], [0.80, 0.35, 0.28],
                [0.59, 0.37, 0.27], [0.80, 0.56, 0.28], [0.84, 0.35, 0.27],
                [0.77, 0.53, 0.27], [0.77, 0.63, 0.28], [0.96, 0.56, 0.26], [0.86, 0.52, 0.26]
            ],
            "triangles": [[3, 4, 5], [6, 7, 8], [9, 10, 11], [12, 13, 15], [13, 14, 15]]
        },
        // greeny 1
        {
            "material": { "ambient": [0.804, 0.675, 0.286], "diffuse": [0.2, 0.3, 0.2], "specular": [0.15, 0.15, 0.15], "n": 4 },
            "vertices": [
                [0.58, 0.71, 0.27], [0.63, 0.56, 0.25], [0.72, 0.59, 0.25], [0.7, 0.74, 0.27],
                [0.56, 0.68, 0.26], [0.54, 0.52, 0.26], [0.63, 0.55, 0.26],
                [0.8, 0.45, 0.26], [0.88, 0.51, 0.26], [0.87, 0.35, 0.26],
                [0.1, 0.85, 0.275], [0.28, 0.14, 0.28], [0.14, 0.84, 0.28], [0.30, 0.14, 0.28],
                [0.12, 0.09, 0.27], [0.28, 0.19, 0.27], [0.44, 0.09, 0.27]
            ],
            "triangles": [[1, 0, 2], [0, 3, 2], [5, 4, 6], [7, 8, 9], [11, 10, 12], [11, 12, 13], [14, 15, 16]]
        },
        // cyan
        {
            "material": { "ambient": [0.263, 0.467, 0.325], "diffuse": [0.2, 0.4, 0.3], "specular": [0.25, 0.25, 0.25], "n": 6 },
            "vertices": [
                [0.58, 0.37, 0.25], [0.82, 0.36, 0.28], [0.81, 0.29, 0.28], [0.62, 0.25, 0.25],
                [0.45, 0.36, 0.22], [0.72, 0.41, 0.20], [0.73, 0.30, 0.20], [0.64, 0.27, 0.20],
                [0.54, 0.13, 0.22], [0.58, 0.37, 0.22], [0.65, 0.37, 0.22], [0.68, 0.12, 0.22],
                [-0.15, 0.73, 0.26], [0.03, 1.2, 0.28], [0.43, 0.91, 0.24], [0.94, 1.09, 0.28]

            ],
            "triangles": [[0, 1, 3], [3, 1, 2], [4, 5, 7], [7, 5, 6], [8, 9, 11], [9, 10, 11], [12, 13, 14], [13, 15, 14]]
        },
        // brown
        {
            "material": { "ambient": [0.275, 0.099, 0.106], "diffuse": [0.0, 0.0, 0.0], "specular": [0, 0, 0], "n": 5 },
            "vertices": [
                [0.60, 0.78, 0.285], [0.665, 0.79, 0.285], [0.49, 0.60, 0.285], [0.68, 0.54, 0.285], [0.82, 0.66, 0.285],
                [0.57, 0.71, 0.26], [0.63, 0.74, 0.26], [0.6, 0.64, 0.26],
                [0.65, 0.74, 0.26], [0.70, 0.75, 0.26], [0.72, 0.66, 0.26],
                [0.60, 0.26, 0.18], [0.73, 0.38, 0.18], [0.74, 0.29, 0.18],
                [0.51, 0.08, 0.23], [0.65, 0.19, 0.23], [0.66, 0.09, 0.23],
                [-0.02, 0.48, 0.25], [0.56, 0.43, 0.25], [0.55, 0.42, 0.25], [-0.04, 0.45, 0.25],
                [0.64, 0.79, 0.27], [0.64, 0.84, 0.27], [0.67, 0.85, 0.27]
            ],
            "triangles": [[0, 2, 4], [1, 2, 4], [2, 3, 4], [5, 6, 7], [8, 9, 10], [11, 12, 13], [14, 15, 16], [17, 18, 19], [17, 20, 19], [21, 22, 23]]
        },
        // white 2
        {
            "material": { "ambient": [0.941, 0.808, 0.412], "diffuse": [0.3, 0.3, 0.3], "specular": [0.25, 0.25, 0.25], "n": 4 },
            "vertices": [
                [0.43, 0.43, 0.26], [0.53, 0.61, 0.25], [0.61, 0.54, 0.25],
                [0.96, 0.56, 0.26], [0.86, 0.52, 0.26], [0.78, 0.45, 0.25], [0.89, 0.42, 0.25],
                [0.608, 0.642, 0.26], [0.619, 0.661, 0.26], [0.642, 0.661, 0.26],
                [0.668, 0.668, 0.26], [0.684, 0.681, 0.26], [0.704, 0.674, 0.26],
                [0.621, 0.61, 0.26], [0.707, 0.604, 0.26], [0.665, 0.584, 0.26]
            ],
            "triangles": [[0, 1, 2], [4, 3, 6], [4, 6, 5], [7, 8, 9], [10, 11, 12], [13, 14, 15]]
        },
        // boarder
        {
            "material": { "ambient": [0.314, 0.149, 0.125], "diffuse": [0.0, 0.0, 0.0], "specular": [0.0, 0.0, 0.0], "n": 5 },
            "vertices": [
                [-0.1, 1.1, 0.06], [1.1, 1.1, 0.06], [0.85, 0.86, 0.06], [0.148, 0.85, 0.06],
                [-0.1, -0.1, 0.06], [0.153, 0.149, 0.06], [0.854, 0.15, 0.06], [1.1, -0.1, 0.06]
            ],
            "triangles": [[0, 1, 2], [0, 2, 3], [4, 5, 6], [4, 6, 7], [4, 5, 0], [0, 3, 5], [1, 2, 6], [6, 1, 7]]
        },
        // boarder 2
        {
            "material": { "ambient": [1, 1, 1], "diffuse": [0.0, 0.0, 0.0], "specular": [0.0, 0.0, 0.0], "n": 5 },
            "vertices": [
                [-0.1, 1.1, 0.01], [1.1, 1.1, 0.01], [0.85, 0.86, 0.01], [0.148, 0.85, 0.01],
                [-0.1, -0.1, 0.01], [0.153, 0.149, 0.01], [0.854, 0.15, 0.01], [1.1, -0.1, 0.01]
            ],
            "triangles": [[0, 1, 2], [0, 2, 3], [4, 5, 6], [4, 6, 7], [4, 5, 0], [0, 3, 5], [1, 2, 6], [6, 1, 7]]
        },
        // greeny 2
        {
            "material": { "ambient": [0.745, 0.580, 0.255], "diffuse": [0.2, 0.3, 0.2], "specular": [0.15, 0.15, 0.15], "n": 5 },
            "vertices": [
                [0.65, 0.62, 0.25], [0.657, 0.667, 0.25], [0.686, 0.63, 0.25],
                [0.65, 0.57, 0.26], [0.71, 0.59, 0.26], [0.69, 0.51, 0.26]

            ],
            "triangles": [[0, 1, 2], [3, 4, 5]]
        }
    ];

    if (showCustom) {
        inputTriangles = inputTriangles_custom;
    }

    

    var w = context.canvas.width;
    var h = context.canvas.height;
    var imagedata = context.createImageData(w, h);

    var eyePos = new Vector(0.5, 0.5, -0.5);
    var viewUpVec = new Vector(0, 1, 0);
    var lookAtVec = new Vector(0, 0, 1);
    var view = { eye: eyePos, at: lookAtVec, up: viewUpVec };
    

    var windowDist = 0.5;
    var windowCenter = new Vector(0.5, 0.5, 0);

    var eyePointSlope;
    var denom, isectT;
    var ctrToIsect;

    var planeCenter = Vector.add(view.eye, Vector.normalize(view.at));
    //console.log("windowCenter: " + windowCenter.toString());
    var planeD = -Vector.dot(view.at, planeCenter);
    var num = -Vector.dot(view.at, view.eye) - planeD;

    var planeX = Vector.normalize(Vector.cross(view.up, view.at));
    var planeY = Vector.normalize(Vector.cross(planeX,view.at));


    var lightPos = new Vector(-3, 1, -0.5);
    var lightColor = [1, 1, 1]; // default: [1, 1, 1]
    var light = { pos: lightPos, color: lightColor};



    if (inputTriangles != String.null) {
        var c = new Color(0, 0, 0, 0); // the color at the pixel: black
       
        var cb = new Color(0, 0, 0, 255); // the color for empty pixels: black
        //var shadowOcc = 1; // 0 = occluded, 1 = not occluded
        //var n = inputTriangles.length;
        //console.log("number of files: " + n);

        var point;

        // goes through each canvas pixel, and only renders those which intersect a tiangle
        for (var y = 0.0; y < h; y++) {
            for (var x = 0.0; x < w; x++) {
                point = new Vector(x / w, y / h, 0);

                var surface = null;
           

                //view stuff
                eyePointSlope = Vector.sub(point, view.eye);
                denom = Vector.dot(view.at, eyePointSlope);

                isectT = num / denom;
                ctrToIsect = Vector.sub(Vector.add(view.eye, Vector.multConst(eyePointSlope, isectT)), planeCenter);
                var px = Vector.dot(planeX, ctrToIsect) * w / 2 + w / 2;
                var py = Vector.dot(planeY, ctrToIsect) * h / 2 + h / 2;
                //var pz = Vector3D.dot(planeY, ctrToIsect) * h / 2 + h / 2;


                // find closest surface intersection
                surface = surfaceRayCast(inputTriangles, point, eyePointSlope, context);
                
                if (surface != null) {

                    var shaded = shade(surface, 0, view.eye, light);

                    c.change(
                        Math.min(shaded[0] * 255, 255),
                        Math.min(shaded[1] * 255, 255),
                        Math.min(shaded[2] * 255, 255),
                        Math.min(shaded[3] * 255, 255));
                    //console.log("c: " + c.toString());

                    drawPixel(imagedata, px, py, c);

                }
                else { // no surface intersection
                    // draw background (black by default)
                    drawPixel(imagedata, px, py, cb);

                }
            }
        }
        console.log("success!");
        context.putImageData(imagedata, 0, 0);
    } // end if triangle files found
} // end draw input triangles


//draw 2d projections traingle from the JSON file at class github
function drawInputTrainglesUsingPaths(context) {
    var inputTriangles = getInputTriangles();
    
    if (inputTriangles != String.null) { 
        var c = new Color(0,0,0,0); // the color at the pixel: black
        var w = context.canvas.width;
        var h = context.canvas.height;
        var n = inputTriangles.length; 
        //console.log("number of files: " + n);

        // Loop over the input files
        for (var f=0; f<n; f++) {
        	var tn = inputTriangles[f].triangles.length;
        	//console.log("number of triangles in this files: " + tn);
        	
        	// Loop over the triangles, draw each in 2d
        	for(var t=0; t<tn; t++){
        		var vertex1 = inputTriangles[f].triangles[t][0];
        		var vertex2 = inputTriangles[f].triangles[t][1];
        		var vertex3 = inputTriangles[f].triangles[t][2];

        		var vertexPos1 = inputTriangles[f].vertices[vertex1];
        		var vertexPos2 = inputTriangles[f].vertices[vertex2];
        		var vertexPos3 = inputTriangles[f].vertices[vertex3];
        		//console.log("vertexPos1 " + vertexPos1);
        		//console.log("vertexPos2 " + vertexPos2);
        		//console.log("vertexPos3 " + vertexPos3);
        		
            	context.fillStyle = 
            	    "rgb(" + Math.floor(inputTriangles[f].material.diffuse[0]*255)
            	    +","+ Math.floor(inputTriangles[f].material.diffuse[1]*255)
            	    +","+ Math.floor(inputTriangles[f].material.diffuse[2]*255) +")"; // diffuse color
            
            	var path=new Path2D();
            	path.moveTo(w*vertexPos1[0],h*vertexPos1[1]);
            	path.lineTo(w*vertexPos2[0],h*vertexPos2[1]);
            	path.lineTo(w*vertexPos3[0],h*vertexPos3[1]);
            	path.closePath();
            	context.fill(path);

        	} // end for triangles
        } // end for files
    } // end if triangle files found
} // end draw input triangles

// put random points in the boxes from the class github
function drawRandPixelsInInputBoxes(context) {
    var inputBoxes = getInputBoxes();
    var w = context.canvas.width;
    var h = context.canvas.height;
    var imagedata = context.createImageData(w,h);
    const PIXEL_DENSITY = 0.1;
    var numCanvasPixels = (w*h)*PIXEL_DENSITY; 
    
    if (inputBoxes != String.null) { 
	    var x  = 0; var y  = 0; // pixel coord init
        var lx = 0; var rx = 0; // input lx, rx from boxes.json
        var by = 0; var ty = 0; // input by, ty from boxes.json
        var fz = 0; var rz = 0; // input fz, rz from boxes.json
        var numBoxPixels = 0; // init num pixels in boxes
        var c = new Color(0,0,0,0); // init the box color
        var n = inputBoxes.length; // the number of input boxes
        //console.log("number of ellipses: " + n);

        // Loop over the ellipsoids, draw rand pixels in each
        for (var b=0; b<n; b++) {
			// input lx,rx,by,ty on canvas
			lx = w*inputBoxes[b].lx;
			rx = w*inputBoxes[b].rx;
			by = h*inputBoxes[b].by;
			ty = h*inputBoxes[b].ty;           
			
            numBoxesPixels  = (rx-lx)*(ty-by); // projected box area 
            numBoxesPixels *= PIXEL_DENSITY;  // percentage of box area to render to pixels
            numBoxesPixels  = Math.round(numBoxesPixels);
           
            //console.log("num box pixels: "+numBoxesPixels);
            
			c.change(
                inputBoxes[b].diffuse[0]*255,
                inputBoxes[b].diffuse[1]*255,
                inputBoxes[b].diffuse[2]*255,
                255); // box diffuse color
            for (var p=0; p<numBoxesPixels; p++) {
                do {
                    x = Math.floor(Math.random()*w); 
                    y = Math.floor(Math.random()*h); 
                } while ( x<lx || x>rx || y>ty || y<by ) // inside the projection
                drawPixel(imagedata,x,y,c);
                //console.log("color: ("+c.r+","+c.g+","+c.b+")");
                //console.log("x: " + x);
                //console.log("y: " + y);
            } // end for pixels in box
        } // end for boxes
        context.putImageData(imagedata, 0, 0);
    } // end if boxes found
} // end draw rand pixels in input boxes

//draw 2d projections boxes from the JSON file at class github
function drawInputBoxesUsingPaths(context) {
    var inputBoxes = getInputBoxes();
    var n = inputBoxes.length; // the number of input boxes
	
    if (inputBoxes != String.null) { 
		var w = context.canvas.width;
        var h = context.canvas.height;
        var c = new Color(0,0,0,0); // the color at the pixel: black
        var x  = 0; var y  = 0; // pixel coord init
        var lx = 0; var rx = 0; // input lx, rx from boxes.json
        var by = 0; var ty = 0; // input by, ty from boxes.json
        var fz = 0; var rz = 0; // input fz, rz from boxes.json
        //console.log("number of files: " + n);

        // Loop over the input files
        for (var b=0; b<n; b++) {
				
			// input lx,rx,by,ty on canvas
			lx = w*inputBoxes[b].lx;
			rx = w*inputBoxes[b].rx;
			by = h*inputBoxes[b].by;
			ty = h*inputBoxes[b].ty; 
        		
            context.fillStyle = 
            	"rgb(" + Math.floor(inputBoxes[b].diffuse[0]*255)
            	+","+ Math.floor(inputBoxes[b].diffuse[1]*255)
            	+","+ Math.floor(inputBoxes[b].diffuse[2]*255) +")"; // diffuse color
            
            var path=new Path2D();
            path.moveTo(lx,ty);
            path.lineTo(lx,by);
            path.lineTo(rx,by);
			path.lineTo(rx,ty);
            path.closePath();
            context.fill(path);

        } // end for files
    } // end if box files found
} // end draw input boxes

/**
 * Toggles the shown graphic between 2 versions when 'Space' is pressed
 * @param {any} evt event data
 */
function toggleDisplay(evt) {
    var canvas = document.getElementById("viewport");
    var context = canvas.getContext("2d");

    // toggles the showCustom boolean when "Space" is pressed
    if (evt.code == "Space") {
        if (showCustom) {
            showCustom = false;
        }
        else {
            showCustom = true;
        }

        drawInputTraingles(context);
    }
}
/* main -- here is where execution begins after window load */
function main() {
    // Get the canvas and context
    var canvas = document.getElementById("viewport"); 
    var context = canvas.getContext("2d");
 
    // Create the image
    //drawRandPixels(context);
      // shows how to draw pixels

    //drawRandPixelsInInputEllipsoids(context);
      // shows how to draw pixels and read input file

    //drawInputEllipsoidsUsingArcs(context);
      // shows how to read input file, but not how to draw pixels

    //drawRandPixelsInInputTriangles(context);
      // shows how to draw pixels and read input file

    //drawInputTrainglesUsingPaths(context);
      // shows how to read input file, but not how to draw pixels

    window.addEventListener("keydown", toggleDisplay, false);

    drawInputTraingles(context);


    //drawRandPixelsInInputBoxes(context);
      // shows how to draw pixels and read input file
    
    //drawInputBoxesUsingPaths(context);
      // shows how to read input file, but not how to draw pixels
}
