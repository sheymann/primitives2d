var primitives2d = {};

;(function(undefined) {
  'use strict';

  /**
   * Return the euclidian distance between two points of a plane with an
   * orthonormal basis.
   *
   * @param  {number} x1  The X coordinate of the first point.
   * @param  {number} y1  The Y coordinate of the first point.
   * @param  {number} x2  The X coordinate of the second point.
   * @param  {number} y2  The Y coordinate of the second point.
   * @return {number}     The euclidian distance.
   */
  primitives2d.getDistance = function(x0, y0, x1, y1) {
    return Math.sqrt(Math.pow(x1 - x0, 2) + Math.pow(y1 - y0, 2));
  };

  /*******************************************************
   *                   LINES & SEGMENTS
   *******************************************************/

  /**
    * Check if a point is on a line segment.
    * http://stackoverflow.com/a/328122
    *
    * @param  {number} x       The X coordinate of the point to check.
    * @param  {number} y       The Y coordinate of the point to check.
    * @param  {number} x1      The X coordinate of the line start point.
    * @param  {number} y1      The Y coordinate of the line start point.
    * @param  {number} x2      The X coordinate of the line end point.
    * @param  {number} y2      The Y coordinate of the line end point.
    * @param  {number} epsilon The precision (consider the line thickness).
    * @return {boolean}        True if point is "close to" the line segment,
    *                          false otherwise.
  */
  primitives2d.isPointOnSegment = function(x, y, x1, y1, x2, y2, epsilon) {
    var crossProduct = Math.abs((y - y1) * (x2 - x1) - (x - x1) * (y2 - y1)),
        d = sigma.utils.getDistance(x1, y1, x2, y2),
        nCrossProduct = crossProduct / d; // normalized cross product

    return (nCrossProduct < epsilon &&
     Math.min(x1, x2) <= x && x <= Math.max(x1, x2) &&
     Math.min(y1, y2) <= y && y <= Math.max(y1, y2));
  };

  /**
   * Get the a point the line segment [A,B] at a specified distance percentage
   * from the start point.
   *
   * @param  {number} aX The x coorinates of the start point.
   * @param  {number} aY The y coorinates of the start point.
   * @param  {number} bX The x coorinates of the end point.
   * @param  {number} bY The y coorinates of the end point.
   * @param  {number} t  The distance percentage from the start point.
   * @return {object}    The (x,y) coordinates of the point.
   */
  primitives2d.getPointOnLineSegment = function(aX, aY, bX, bY, t) {
    return {
      x: aX + (bX - aX) * t,
      y: aY + (bY - aY) * t
    };
  };

  /**
   * Find the closest point on a line.
   * http://stackoverflow.com/a/328122 & http://stackoverflow.com/a/3120357
   *
   * @param  {number} x  The X coordinate of the point to check.
   * @param  {number} y  The Y coordinate of the point to check.
   * @param  {number} x1 The X coordinate of the line start point.
   * @param  {number} y1 The Y coordinate of the line start point.
   * @param  {number} x2 The X coordinate of the line end point.
   * @param  {number} y2 The Y coordinate of the line end point.
   * @return {x,y}       The coordinates of the closest point.
   */
  primitives2d.getClosestPointToLine = function(x, y, x1, y1, x2, y2) {
    var crossProduct = Math.abs((y - y1) * (x2 - x1) - (x - x1) * (y2 - y1)),
      d = getDistance(x1, y1, x2, y2),
      nCrossProduct = crossProduct / d; // normalized cross product

    // Add the distance to A, moving towards B
    return {
      x: x1 + (x2 - x1) * nCrossProduct,
      y: y1 + (y2 - y1) * nCrossProduct
    }
  };

  /**
   * Find the intersection between two lines, two segments, or one line and one segment.
   * http://jsfiddle.net/justin_c_rounds/Gd2S2/
   *
   * @param  {number} line1x1  The X coordinate of the start point of the first line.
   * @param  {number} line1y1  The Y coordinate of the start point of the first line.
   * @param  {number} line1x2  The X coordinate of the end point of the first line.
   * @param  {number} line1y2  The Y coordinate of the end point of the first line.v
   * @param  {number} line2x1  The X coordinate of the start point of the second line.
   * @param  {number} line2y1  The Y coordinate of the start point of the second line.
   * @param  {number} line2x2  The X coordinate of the end point of the second line.
   * @param  {number} line2y2  The Y coordinate of the end point of the second line.
   * @return {object}           The coordinates of the intersection point.
   */
  primitives2d.getLinesIntersection =
    function(line1x1, line1y1, line1x2, line1y2, line2x1, line2y1, line2x2, line2y2) {
    // if the lines intersect, the result contains the x and y of the intersection
    // (treating the lines as infinite) and booleans for whether line segment 1 or
    // line segment 2 contain the point
    var
      denominator,
      a,
      b,
      numerator1,
      numerator2,
      result = {
        x: null,
        y: null,
        onLine1: false,
        onLine2: false
    };

    denominator =
      ((line2y2 - line2y1) * (line1x2 - line1x1)) -
      ((line2x2 - line2x1) * (line1y2 - line1y1));

    if (denominator == 0) {
        return result;
    }

    a = line1y1 - line2y1;
    b = line1x1 - line2x1;

    numerator1 = ((line2x2 - line2x1) * a) - ((line2y2 - line2y1) * b);
    numerator2 = ((line1x2 - line1x1) * a) - ((line1y2 - line1y1) * b);

    a = numerator1 / denominator;
    b = numerator2 / denominator;

    // if we cast these lines infinitely in both directions, they intersect here:
    result.x = line1x1 + (a * (line1x2 - line1x1));
    result.y = line1y1 + (a * (line1y2 - line1y1));
    /*
    // it is worth noting that this should be the same as:
      x = line2x1 + (b * (line2x2 - line2x1));
      y = line2x1 + (b * (line2y2 - line2y1));
    */
    // if line1 is a segment and line2 is infinite, they intersect if:
    if (a > 0 && a < 1) {
        result.onLine1 = true;
    }
    // if line2 is a segment and line1 is infinite, they intersect if:
    if (b > 0 && b < 1) {
        result.onLine2 = true;
    }
    // if line1 and line2 are segments, they intersect if both of the above are true
    return result;
  };


  /*******************************************************
   *                      CIRCLES
   *******************************************************/

  /**
   * Return the coordinates of the intersection points of two circles.
   * http://stackoverflow.com/a/12219802
   *
   * @param  {number} x0  The X coordinate of center location of the first
   *                      circle.
   * @param  {number} y0  The Y coordinate of center location of the first
   *                      circle.
   * @param  {number} r0  The radius of the first circle.
   * @param  {number} x1  The X coordinate of center location of the second
   *                      circle.
   * @param  {number} y1  The Y coordinate of center location of the second
   *                      circle.
   * @param  {number} r1  The radius of the second circle.
   * @return {object}     The (xi,yi) coordinates of the intersection points.
   */
  primitives2d.getCircleIntersection = function(x0, y0, r0, x1, y1, r1) {
    var a, dx, dy, d, h, rx, ry, x2, y2;

    // dx and dy are the vertical and horizontal distances between the circle
    // centers:
    dx = x1 - x0;
    dy = y1 - y0;

    // Determine the straight-line distance between the centers:
    d = Math.sqrt((dy * dy) + (dx * dx));

    // Check for solvability:
    if (d > (r0 + r1)) {
        // No solution. circles do not intersect.
        return false;
    }
    if (d < Math.abs(r0 - r1)) {
        // No solution. one circle is contained in the other.
        return false;
    }

    //'point 2' is the point where the line through the circle intersection
    // points crosses the line between the circle centers.

    // Determine the distance from point 0 to point 2:
    a = ((r0 * r0) - (r1 * r1) + (d * d)) / (2.0 * d);

    // Determine the coordinates of point 2:
    x2 = x0 + (dx * a / d);
    y2 = y0 + (dy * a / d);

    // Determine the distance from point 2 to either of the intersection
    // points:
    h = Math.sqrt((r0 * r0) - (a * a));

    // Determine the offsets of the intersection points from point 2:
    rx = -dy * (h / d);
    ry = dx * (h / d);

    // Determine the absolute intersection points:
    var xi = x2 + rx;
    var xi_prime = x2 - rx;
    var yi = y2 + ry;
    var yi_prime = y2 - ry;

    return {xi: xi, xi_prime: xi_prime, yi: yi, yi_prime: yi_prime};
  };

  /*******************************************************
   *                      VECTORS
   *******************************************************/

  /**
   * Get the angle of the vector (in radian).
   *
   * @param  {object} v  The 2d vector with x,y coordinates.
   * @return {number}    The angle of the vector  (in radian).
   */
  primitives2d.getVectorAngle = function(v) {
    return Math.acos( v.x / Math.sqrt(v.x * v.x + v.y * v.y) );
  };

  /**
   * Get the normal vector of the line segment, i.e. the vector
   * orthogonal to the line.
   * http://stackoverflow.com/a/1243614/
   *
   * @param  {number} aX The x coorinates of the start point.
   * @param  {number} aY The y coorinates of the start point.
   * @param  {number} bX The x coorinates of the end point.
   * @param  {number} bY The y coorinates of the end point.
   * @return {object}    The 2d vector with (xi,yi), (xi_prime,yi_prime) coordinates.
   */
  primitives2d.getNormalVector = function(aX, aY, bX, bY) {
    return {
      xi:       -(bY - aY),
      yi:         bX - aX,
      xi_prime:   bY - aY,
      yi_prime: -(bX - aX)
    };
  };

  /**
   * Get the normalized vector.
   *
   * @param  {object} v      The 2d vector with (xi,yi), (xi_prime,yi_prime) coordinates.
   * @param  {number} length The vector length.
   * @return {object}        The normalized vector
   */
  primitives2d.getNormalizedVector = function(v, length) {
    return {
      x: (v.xi_prime - v.xi) / length,
      y: (v.yi_prime - v.yi) / length,
    };
  };

  /*******************************************************
   *              QUADRATIC BEZIER CURVES
   *******************************************************/

  /**
   * Return the control point coordinates for a quadratic bezier curve.
   *
   * @param  {number} x1  The X coordinate of the start point.
   * @param  {number} y1  The Y coordinate of the start point.
   * @param  {number} x2  The X coordinate of the end point.
   * @param  {number} y2  The Y coordinate of the end point.
   * @return {object}     The (x,y) control point coordinates.
   */
  primitives2d.getQuadraticControlPoint = function(x1, y1, x2, y2) {
    return {
      x: (x1 + x2) / 2 + (y2 - y1) / 4,
      y: (y1 + y2) / 2 + (x1 - x2) / 4
    };
  };

  /**
    * Compute the coordinates of the point positioned
    * at length t in the quadratic bezier curve.
    * http://stackoverflow.com/a/5634528
    *
    * @param  {number} t  In [0,1] the step percentage to reach
    *                     the point in the curve from the context point.
    * @param  {number} x1 The X coordinate of the context point.
    * @param  {number} y1 The Y coordinate of the context point.
    * @param  {number} x2 The X coordinate of the ending point.
    * @param  {number} y2 The Y coordinate of the ending point.
    * @param  {number} xi The X coordinate of the control point.
    * @param  {number} yi The Y coordinate of the control point.
    * @return {object}    The (x,y) coordinates of the point.
  */
  primitives2d.getPointOnQuadraticCurve = function(t, x1, y1, x2, y2, xi, yi) {
    return {
      x: Math.pow(1 - t, 2) * x1 + 2 * (1 - t) * t * xi + Math.pow(t, 2) * x2,
      y: Math.pow(1 - t, 2) * y1 + 2 * (1 - t) * t * yi + Math.pow(t, 2) * y2
    };
  };

  /**
    * Check if a point is on a quadratic bezier curve segment with a thickness.
    *
    * @param  {number} x       The X coordinate of the point to check.
    * @param  {number} y       The Y coordinate of the point to check.
    * @param  {number} x1      The X coordinate of the curve start point.
    * @param  {number} y1      The Y coordinate of the curve start point.
    * @param  {number} x2      The X coordinate of the curve end point.
    * @param  {number} y2      The Y coordinate of the curve end point.
    * @param  {number} cpx     The X coordinate of the curve control point.
    * @param  {number} cpy     The Y coordinate of the curve control point.
    * @param  {number} epsilon The precision (consider the line thickness).
    * @return {boolean}        True if (x,y) is on the curve segment,
    *                          false otherwise.
  */
  primitives2d.isPointOnQuadraticCurve = function(x, y, x1, y1, x2, y2, cpx, cpy, epsilon) {
    // Fails if the point is too far from the extremities of the segment,
    // preventing for more costly computation:
    var dP1P2 = sigma.utils.getDistance(x1, y1, x2, y2);
    if (Math.abs(x - x1) > dP1P2 || Math.abs(y - y1) > dP1P2) {
      return false;
    }

    var dP1 = sigma.utils.getDistance(x, y, x1, y1),
        dP2 = sigma.utils.getDistance(x, y, x2, y2),
        t = 0.5,
        r = (dP1 < dP2) ? -0.01 : 0.01,
        rThreshold = 0.001,
        i = 100,
        pt = sigma.utils.getPointOnQuadraticCurve(t, x1, y1, x2, y2, cpx, cpy),
        dt = sigma.utils.getDistance(x, y, pt.x, pt.y),
        old_dt;

    // This algorithm minimizes the distance from the point to the curve. It
    // find the optimal t value where t=0 is the start point and t=1 is the end
    // point of the curve, starting from t=0.5.
    // It terminates because it runs a maximum of i interations.
    while (i-- > 0 &&
      t >= 0 && t <= 1 &&
      (dt > epsilon) &&
      (r > rThreshold || r < -rThreshold)) {
      old_dt = dt;
      pt = sigma.utils.getPointOnQuadraticCurve(t, x1, y1, x2, y2, cpx, cpy);
      dt = sigma.utils.getDistance(x, y, pt.x, pt.y);

      if (dt > old_dt) {
        // not the right direction:
        // halfstep in the opposite direction
        r = -r / 2;
        t += r;
      }
      else if (t + r < 0 || t + r > 1) {
        // oops, we've gone too far:
        // revert with a halfstep
        r = r / 2;
        dt = old_dt;
      }
      else {
        // progress:
        t += r;
      }
    }

    return dt < epsilon;
  };

  /*******************************************************
   *                CUBIC BEZIER CURVES
   *******************************************************/

  /**
    * Compute the coordinates of the point positioned at length t in the cubic
    * bezier curve.
    * http://stackoverflow.com/a/15397596
    *
    * @param  {number} t  In [0,1] the step percentage to reach the point in
    *                     the curve from the context point.
    * @param  {number} x1 The X coordinate of the context point.
    * @param  {number} y1 The Y coordinate of the context point.
    * @param  {number} x2 The X coordinate of the end point.
    * @param  {number} y2 The Y coordinate of the end point.
    * @param  {number} cx The X coordinate of the first control point.
    * @param  {number} cy The Y coordinate of the first control point.
    * @param  {number} dx The X coordinate of the second control point.
    * @param  {number} dy The Y coordinate of the second control point.
    * @return {object}    The (x,y) coordinates of the point at t.
  */
  primitives2d.getPointOnBezierCurve = function(t, x1, y1, x2, y2, cx, cy, dx, dy) {
    // Blending functions:
    var B0_t = Math.pow(1 - t, 3),
        B1_t = 3 * t * Math.pow(1 - t, 2),
        B2_t = 3 * Math.pow(t, 2) * (1 - t),
        B3_t = Math.pow(t, 3);

    return {
      x: (B0_t * x1) + (B1_t * cx) + (B2_t * dx) + (B3_t * x2),
      y: (B0_t * y1) + (B1_t * cy) + (B2_t * dy) + (B3_t * y2)
    };
  };

  /**
    * Check if a point is on a cubic bezier curve segment with a thickness.
    *
    * @param  {number} x       The X coordinate of the point to check.
    * @param  {number} y       The Y coordinate of the point to check.
    * @param  {number} x1      The X coordinate of the curve start point.
    * @param  {number} y1      The Y coordinate of the curve start point.
    * @param  {number} x2      The X coordinate of the curve end point.
    * @param  {number} y2      The Y coordinate of the curve end point.
    * @param  {number} cpx1    The X coordinate of the 1st curve control point.
    * @param  {number} cpy1    The Y coordinate of the 1st curve control point.
    * @param  {number} cpx2    The X coordinate of the 2nd curve control point.
    * @param  {number} cpy2    The Y coordinate of the 2nd curve control point.
    * @param  {number} epsilon The precision (consider the line thickness).
    * @return {boolean}        True if (x,y) is on the curve segment, false otherwise.
  */
  primitives2d.isPointOnBezierCurve = function(x, y, x1, y1, x2, y2, cpx1, cpy1, cpx2, cpy2, epsilon) {
    // Fails if the point is too far from the extremities of the segment,
    // preventing for more costly computation:
    var dP1CP1 = sigma.utils.getDistance(x1, y1, cpx1, cpy1);
    if (Math.abs(x - x1) > dP1CP1 || Math.abs(y - y1) > dP1CP1) {
      return false;
    }

    var dP1 = sigma.utils.getDistance(x, y, x1, y1),
        dP2 = sigma.utils.getDistance(x, y, x2, y2),
        t = 0.5,
        r = (dP1 < dP2) ? -0.01 : 0.01,
        rThreshold = 0.001,
        i = 100,
        pt = sigma.utils.getPointOnBezierCurve(
          t, x1, y1, x2, y2, cpx1, cpy1, cpx2, cpy2),
        dt = sigma.utils.getDistance(x, y, pt.x, pt.y),
        old_dt;

    // This algorithm minimizes the distance from the point to the curve. It
    // find the optimal t value where t=0 is the start point and t=1 is the end
    // point of the curve, starting from t=0.5.
    // It terminates because it runs a maximum of i interations.
    while (i-- > 0 &&
      t >= 0 && t <= 1 &&
      (dt > epsilon) &&
      (r > rThreshold || r < -rThreshold)) {
      old_dt = dt;
      pt = sigma.utils.getPointOnBezierCurve(
        t, x1, y1, x2, y2, cpx1, cpy1, cpx2, cpy2);
      dt = sigma.utils.getDistance(x, y, pt.x, pt.y);

      if (dt > old_dt) {
        // not the right direction:
        // halfstep in the opposite direction
        r = -r / 2;
        t += r;
      }
      else if (t + r < 0 || t + r > 1) {
        // oops, we've gone too far:
        // revert with a halfstep
        r = r / 2;
        dt = old_dt;
      }
      else {
        // progress:
        t += r;
      }
    }

    return dt < epsilon;
  };

  /*******************************************************
   *                      UTILS
   *******************************************************/

  /**
   * Scale a value from the range [baseMin, baseMax] to the range [limitMin, limitMax].
   *
   * @param  {number} value    The value to rescale.
   * @param  {number} baseMin  The min value of the range of origin.
   * @param  {number} baseMax  The max value of the range of origin.
   * @param  {number} limitMin The min value of the range of destination.
   * @param  {number} limitMax The max value of the range of destination.
   * @return {number}          The scaled value.
   */
  primitives2d.scaleRange = function(value, baseMin, baseMax, limitMin, limitMax) {
    return ((limitMax - limitMin) * (value - baseMin) / (baseMax - baseMin)) + limitMin;
  };


}).call(this);
