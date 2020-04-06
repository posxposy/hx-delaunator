package delaunator;

@:keep
#if js
@:expose
#elseif (java || cs)
@:nativeGen
#elseif cpp
@:unreflective
#end
final class Delaunator {
	static final EPSILON:Float = Math.pow(2, -52);
	static final EDGE_STACK:IntContainer = new IntContainer(512);

	public static function from(points:Array<Array<Float>>):Delaunator {
		final n = points.length;
		final coords = new SingleContainer(n * 2);
		for (i in 0...points.length) {
			final point = points[i];
			coords[2 * i + 0] = point[0];
			coords[2 * i + 1] = point[1];
		}

		return new Delaunator(coords);
	}

	public static function fromPoints<T:{x:Float, y:Float}>(points:Array<T>):Delaunator {
		final n = points.length;
		final coords = new SingleContainer(n * 2);
		for (i in 0...points.length) {
			final point = points[i];
			coords[2 * i + 0] = point.x;
			coords[2 * i + 1] = point.y;
		}

		return new Delaunator(coords);
	}

	public var triangles(default, null):UIntContainer;
	public var halfedges(default, null):IntContainer;

	public final coords:SingleContainer;

	final hashSize:Int;

	final hullPrev:UIntContainer;
	final hullNext:UIntContainer;
	final hullTri:UIntContainer;
	final ids:UIntContainer;
	final hullHash:IntContainer;
	final dists:FloatContainer;

	final center:Point;
	var hull:IntContainer;
	var trianglesLen:Int;
	var hullStart:Int;

	function new(coords:SingleContainer) {
		final n = coords.length >> 1;

		this.coords = coords;

		// arrays that will store the triangulation graph
		final maxTriangles = Std.int(Math.max(2 * n - 5, 0));
		triangles = new UIntContainer(maxTriangles * 3);
		halfedges = new IntContainer(maxTriangles * 3);

		// temporary arrays for tracking the edges of the advancing convex hull
		hashSize = Math.ceil(Math.sqrt(n));
		hullPrev = new UIntContainer(n); // edge to prev edge
		hullNext = new UIntContainer(n); // edge to next edge
		hullTri = new UIntContainer(n); // edge to adjacent triangle
		hullHash = new IntContainer(hashSize).fill(-1); // angular edge

		// temporary arrays for sorting points
		ids = new UIntContainer(n);
		dists = new FloatContainer(n);

		trianglesLen = 0;

		center = new Point(0, 0);

		update();
	}

	function update():Void {
		final n = coords.length >> 1;

		// populate an array of point indices; calculate input data bbox
		var minX = Math.POSITIVE_INFINITY, minY = Math.POSITIVE_INFINITY;
		var maxX = Math.NEGATIVE_INFINITY, maxY = Math.NEGATIVE_INFINITY;

		for (i in 0...n) {
			final x = coords[2 * i + 0];
			final y = coords[2 * i + 1];
			if (x < minX) {
				minX = x;
			}
			if (y < minY) {
				minY = y;
			}
			if (x > maxX) {
				maxX = x;
			}
			if (y > maxY) {
				maxY = y;
			}

			this.ids[i] = i;
		}

		final cx = (minX + maxX) * 0.5;
		final cy = (minY + maxY) * 0.5;

		var minDist = Math.POSITIVE_INFINITY;
		var i0 = 0, i1 = 0, i2 = 0;

		// pick a seed point close to the center
		for (i in 0...n) {
			final d = dist(cx, cy, coords[2 * i], coords[2 * i + 1]);
			if (d < minDist) {
				i0 = i;
				minDist = d;
			}
		}
		final i0x = coords[2 * i0];
		final i0y = coords[2 * i0 + 1];

		minDist = Math.POSITIVE_INFINITY;

		// find the point closest to the seed
		for (i in 0...n) {
			if (i == i0) {
				continue;
			}
			final d = dist(i0x, i0y, coords[2 * i], coords[2 * i + 1]);
			if (d < minDist && d > 0) {
				i1 = i;
				minDist = d;
			}
		}
		var i1x = coords[2 * i1];
		var i1y = coords[2 * i1 + 1];

		var minRadius = Math.POSITIVE_INFINITY;

		// find the third point which forms the smallest circumcircle with the first two
		for (i in 0...n) {
			if (i == i0 || i == i1) {
				continue;
			}
			final r = circumradius(i0x, i0y, i1x, i1y, coords[2 * i], coords[2 * i + 1]);
			if (r < minRadius) {
				i2 = i;
				minRadius = r;
			}
		}
		var i2x = coords[2 * i2];
		var i2y = coords[2 * i2 + 1];

		if (minRadius == Math.POSITIVE_INFINITY) {
			return;
			// order collinear points by dx (or dy if all x are identical)
			// and return the list as a hull
			for (i in 0...n) {
				final a = coords[2 * i] - coords[0];
				final b = coords[2 * i + 1] - coords[1];
				dists[i] = a != 0 ? a : b;
			}
			quicksort(0, n - 1);
			final hull = new IntContainer(n);
			var j = 0;
			var d0 = Math.NEGATIVE_INFINITY;
			for (i in 0...n) {
				final id = ids[i];
				if (dists[id] > d0) {
					hull[j++] = id;
					d0 = dists[id];
				}
			}

			this.hull = hull.subarray(0, j);
			triangles = new UIntContainer(0);
			halfedges = new IntContainer(0);
			return;
		}

		// swap the order of the seed points for counter-clockwise orientation
		if (orient(i0x, i0y, i1x, i1y, i2x, i2y)) {
			final i = i1;
			final x = i1x;
			final y = i1y;
			i1 = i2;
			i1x = i2x;
			i1y = i2y;
			i2 = i;
			i2x = x;
			i2y = y;
		}

		circumcenter(center, i0x, i0y, i1x, i1y, i2x, i2y);

		for (i in 0...n) {
			dists[i] = dist(coords[2 * i], coords[2 * i + 1], center.x, center.y);
		}

		// sort the points by distance from the seed triangle circumcenter
		quicksort(0, n - 1);

		// set up the seed triangle as the starting hull
		hullStart = i0;
		var hullSize = 3;

		hullNext[i0] = hullPrev[i2] = i1;
		hullNext[i1] = hullPrev[i0] = i2;
		hullNext[i2] = hullPrev[i1] = i0;

		hullTri[i0] = 0;
		hullTri[i1] = 1;
		hullTri[i2] = 2;

		hullHash.fill(-1);
		hullHash[hashKey(i0x, i0y)] = i0;
		hullHash[hashKey(i1x, i1y)] = i1;
		hullHash[hashKey(i2x, i2y)] = i2;

		trianglesLen = 0;
		addTriangle(i0, i1, i2, -1, -1, -1);

		var xp = 0.0;
		var yp = 0.0;

		for (k in 0...ids.length) {
			final i = ids[k];
			final x = coords[2 * i];
			final y = coords[2 * i + 1];

			// skip near-duplicate points
			if (k > 0 && Math.abs(x - xp) <= EPSILON && Math.abs(y - yp) <= EPSILON) {
				continue;
			}
			xp = x;
			yp = y;

			// skip seed triangle points
			if (i == i0 || i == i1 || i == i2) {
				continue;
			}

			// find a visible edge on the convex hull using edge hash
			var start = 0;
			final key = hashKey(x, y);
			for (j in 0...hashSize) {
				start = hullHash[(key + j) % hashSize];
				if (start != -1 && start != hullNext[start]) {
					break;
				}
			}

			start = hullPrev[start];
			var e = start;
			var q = hullNext[e];
			while (!orient(x, y, coords[2 * e], coords[2 * e + 1], coords[2 * q], coords[2 * q + 1])) {
				e = q;
				if (e == start) {
					e = -1;
					break;
				}
				q = hullNext[e];
			}

			if (e == -1) {
				continue; // likely a near-duplicate point; skip it
			}

			// add the first triangle from the point
			var t:Int = addTriangle(e, i, hullNext[e], -1, -1, hullTri[e]);

			// recursively flip triangles from the point until they satisfy the Delaunay condition
			hullTri[i] = legalize(t + 2);
			hullTri[e] = t; // keep track of boundary triangles on the hull
			hullSize++;

			// walk forward through the hull, adding more triangles and flipping recursively
			var next = hullNext[e];
			q = hullNext[next];
			while (orient(x, y, coords[2 * next], coords[2 * next + 1], coords[2 * q], coords[2 * q + 1])) {
				t = addTriangle(next, i, q, hullTri[i], -1, hullTri[next]);
				hullTri[i] = legalize(t + 2);
				hullNext[next] = next; // mark as removed
				hullSize--;
				next = q;
				q = hullNext[next];
			}

			// walk backward from the other side, adding more triangles and flipping
			if (e == start) {
				q = hullPrev[e];
				while (orient(x, y, coords[2 * q], coords[2 * q + 1], coords[2 * e], coords[2 * e + 1])) {
					t = addTriangle(q, i, e, -1, hullTri[e], hullTri[q]);
					legalize(t + 2);
					hullTri[q] = t;
					hullNext[e] = e; // mark as removed
					hullSize--;
					e = q;
					q = hullPrev[e];
				}
			}
			// update the hull indices
			hullStart = hullPrev[i] = e;
			hullNext[e] = hullPrev[next] = i;
			hullNext[i] = next;

			// save the two new edges in the hash table
			hullHash[hashKey(x, y)] = i;
			hullHash[hashKey(coords[2 * e], coords[2 * e + 1])] = e;
		}

		hull = new IntContainer(hullSize);
		var e = hullStart;
		for (i in 0...hullSize) {
			hull[i] = e;
			e = hullNext[e];
		}
		// trim typed triangle mesh arrays
		triangles = triangles.subarray(0, trianglesLen);
		halfedges = halfedges.subarray(0, trianglesLen);
	}

	function legalize(a:Int):Int {
		var i = 0;
		var ar = 0;

		// recursion eliminated with a fixed-size stack
		while (true) {
			final b = halfedges[a];
			/* if the pair of triangles doesn't satisfy the Delaunay condition
			 * (p1 is inside the circumcircle of [p0, pl, pr]), flip them,
			 * then do the same check/flip recursively for the new pair of triangles
			 *
			 *           pl                    pl
			 *          /||\                  /  \
			 *       al/ || \bl            al/    \a
			 *        /  ||  \              /      \
			 *       /  a||b  \    flip    /___ar___\
			 *     p0\   ||   /p1   =>   p0\---bl---/p1
			 *        \  ||  /              \      /
			 *       ar\ || /br             b\    /br
			 *          \||/                  \  /
			 *           pr                    pr
			 */
			final a0 = a - a % 3;
			ar = a0 + (a + 2) % 3;

			if (b == -1) { // convex hull edge
				if (i == 0) {
					break;
				}
				a = EDGE_STACK[--i];
				continue;
			}

			final b0 = b - b % 3;
			final al = a0 + (a + 1) % 3;
			final bl = b0 + (b + 2) % 3;

			final p0 = triangles[ar];
			final pr = triangles[a];
			final pl = triangles[al];
			final p1 = triangles[bl];

			final illegal = inCircle(coords[2 * p0], coords[2 * p0 + 1], coords[2 * pr], coords[2 * pr + 1], coords[2 * pl], coords[2 * pl + 1],
				coords[2 * p1], coords[2 * p1 + 1]);

			if (illegal) {
				triangles[a] = p1;
				triangles[b] = p0;

				final hbl = halfedges[bl];

				// edge swapped on the other side of the hull (rare); fix the halfedge reference
				if (hbl == -1) {
					var e = hullStart;
					do {
						if (hullTri[e] == bl) {
							hullTri[e] = a;
							break;
						}
						e = hullPrev[e];
					} while (e != hullStart);
				}

				link(a, hbl);
				link(b, halfedges[ar]);
				link(ar, bl);

				final br = b0 + (b + 1) % 3;

				// don't worry about hitting the cap: it can only happen on extremely degenerate input
				if (i < EDGE_STACK.length) {
					EDGE_STACK[i++] = br;
				}
			} else {
				if (i == 0) {
					break;
				}
				a = EDGE_STACK[--i];
			}
		}

		return ar;
	}

	function quicksort(left:Int, right:Int):Void {
		if (right - left <= 20) {
			for (i in (left + 1)...(right + 1)) {
				final temp = ids[i];
				final tempDist = dists[temp];
				var j = i - 1;
				while (j >= left && dists[ids[j]] > tempDist) {
					ids[j + 1] = ids[j--];
				}
				ids[j + 1] = temp;
			}
		} else {
			final median = (left + right) >> 1;
			var i = left + 1;
			var j = right;
			swap(ids, median, i);
			if (dists[ids[left]] > dists[ids[right]]) {
				swap(ids, left, right);
			}
			if (dists[ids[i]] > dists[ids[right]]) {
				swap(ids, i, right);
			}
			if (dists[ids[left]] > dists[ids[i]]) {
				swap(ids, left, i);
			}

			final temp = ids[i];
			final tempDist = dists[temp];
			while (true) {
				do {
					i++;
				} while (dists[ids[i]] < tempDist);
				do {
					j--;
				} while (dists[ids[j]] > tempDist);
				if (j < i) {
					break;
				}
				swap(ids, i, j);
			}
			ids[left + 1] = ids[j];
			ids[j] = temp;

			if (right - i + 1 >= j - left) {
				quicksort(i, right);
				quicksort(left, j - 1);
			} else {
				quicksort(left, j - 1);
				quicksort(i, right);
			}
		}
	}

	function circumradius(ax:Float, ay:Float, bx:Float, by:Float, cx:Float, cy:Float):Float {
		final dx = bx - ax;
		final dy = by - ay;
		final ex = cx - ax;
		final ey = cy - ay;

		final bl = dx * dx + dy * dy;
		final cl = ex * ex + ey * ey;
		final d = 0.5 / (dx * ey - dy * ex);

		final x = (ey * bl - dy * cl) * d;
		final y = (dx * cl - ex * bl) * d;

		return x * x + y * y;
	}

	function circumcenter(point:Point, ax:Float, ay:Float, bx:Float, by:Float, cx:Float, cy:Float):Void {
		final dx = bx - ax;
		final dy = by - ay;
		final ex = cx - ax;
		final ey = cy - ay;

		final bl = dx * dx + dy * dy;
		final cl = ex * ex + ey * ey;
		final d = 0.5 / (dx * ey - dy * ex);

		point.x = ax + (ey * bl - dy * cl) * d;
		point.y = ay + (dx * cl - ex * bl) * d;
	}

	inline function inCircle(ax:Float, ay:Float, bx:Float, by:Float, cx:Float, cy:Float, px:Float, py:Float):Bool {
		final dx = ax - px;
		final dy = ay - py;
		final ex = bx - px;
		final ey = by - py;
		final fx = cx - px;
		final fy = cy - py;

		final ap = dx * dx + dy * dy;
		final bp = ex * ex + ey * ey;
		final cp = fx * fx + fy * fy;

		return dx * (ey * cp - bp * fy) - dy * (ex * cp - bp * fx) + ap * (ex * fy - ey * fx) < 0;
	}

	// monotonically increases with real angle, but doesn't need expensive trigonometry
	inline function pseudoAngle(dx:Float, dy:Float) {
		final p = dx / (Math.abs(dx) + Math.abs(dy));
		return (dy > 0 ? 3 - p : 1 + p) / 4; // [0..1]
	}

	inline function hashKey(x:Float, y:Float) {
		return Math.floor(pseudoAngle(x - center.x, y - center.y) * hashSize) % hashSize;
	}

	inline function link(a:Int, b:Int):Void {
		halfedges[a] = b;
		if (b != -1) {
			halfedges[b] = a;
		}
	}

	// add a new triangle given vertex indices and adjacent half-edge ids
	inline function addTriangle(i0:Int, i1:Int, i2:Int, a:Int, b:Int, c:Int):Int {
		final t = trianglesLen;

		triangles[t] = i0;
		triangles[t + 1] = i1;
		triangles[t + 2] = i2;

		link(t, a);
		link(t + 1, b);
		link(t + 2, c);

		trianglesLen += 3;

		return t;
	}

	// return 2d orientation sign if we're confident in it through J. Shewchuk's error bound check
	inline function orientIfSure(px:Float, py:Float, rx:Float, ry:Float, qx:Float, qy:Float):Float {
		final l = (ry - py) * (qx - px);
		final r = (rx - px) * (qy - py);
		return Math.abs(l - r) >= 3.3306690738754716e-16 * Math.abs(l + r) ? l - r : 0;
	}

	// a more robust orientation test that's stable in a given triangle (to fix robustness issues)
	inline function orient(rx:Float, ry:Float, qx:Float, qy:Float, px:Float, py:Float):Bool {
		inline function check(a:Float, b:Float):Float {
			return a != 0 ? a : b;
		}
		return check(check(orientIfSure(px, py, rx, ry, qx, qy), orientIfSure(rx, ry, qx, qy, px, py)), orientIfSure(qx, qy, px, py, rx, ry)) < 0;
	}

	inline function swap(arr:UIntContainer, i:Int, j:Int):Void {
		final tmp = arr[i];
		arr[i] = arr[j];
		arr[j] = tmp;
	}

	inline function dist(ax:Float, ay:Float, bx:Float, by:Float):Float {
		final dx = ax - bx;
		final dy = ay - by;
		return dx * dx + dy * dy;
	}
}

private final class Point {
	public var x:Float;
	public var y:Float;

	public inline function new(x:Float, y:Float) {
		this.x = x;
		this.y = y;
	}
}

private typedef RealFloat = #if sys Single #else Float #end;
private typedef RealUInt = #if sys UInt #else Int #end;
private typedef NativeContainer<T> = haxe.ds.Vector<T>; // todo: std::vector for cpp
private typedef UIntContainerBase = #if js js.lib.Uint32Array #else NativeContainer<RealUInt> #end;
private typedef IntContainerBase = #if js js.lib.Int32Array #else NativeContainer<Int> #end;
private typedef SingleContainerBase = #if js js.lib.Float32Array #else NativeContainer<RealFloat> #end;
private typedef FloatContainerBase = #if js js.lib.Float64Array #else NativeContainer<Float> #end;

@:forward(length)
private abstract UIntContainer(UIntContainerBase) from UIntContainerBase {
	public inline function new(size:Int) {
		#if js
		this = new js.lib.Uint32Array(size);
		#else
		this = new NativeContainer<RealUInt>(size);
		#end
	}

	@:arrayAccess
	public inline function get(i:Int):RealUInt {
		return this[i];
	}

	@:arrayAccess
	public inline function set(i:Int, v:RealUInt):RealUInt {
		this[i] = v;
		return v;
	}

	public inline function subarray(from:Int, to:Int):UIntContainer {
		#if js
		return this.subarray(from, to);
		#else
		return NativeContainerTools.subarray(this, from, to);
		#end
	}
}

@:forward(length)
private abstract IntContainer(IntContainerBase) from IntContainerBase {
	public inline function new(size:Int) {
		#if js
		this = new js.lib.Int32Array(size);
		#else
		this = new NativeContainer<Int>(size);
		#end
	}

	@:arrayAccess
	public inline function get(i:Int):Int {
		return this[i];
	}

	@:arrayAccess
	public inline function set(i:Int, v:Int):Int {
		this[i] = v;
		return v;
	}

	public inline function subarray(from:Int, to:Int):IntContainer {
		#if js
		return this.subarray(from, to);
		#else
		return NativeContainerTools.subarray(this, from, to);
		#end
	}

	public inline function fill(value:Int):IntContainer {
		#if js
		return this.fill(value);
		#else
		return NativeContainerTools.fill(this, value);
		#end
	}
}

@:forward(length)
private abstract SingleContainer(SingleContainerBase) {
	public inline function new(size:Int) {
		#if js
		this = new js.lib.Float32Array(size);
		#elseif sys
		this = new NativeContainer(size);
		#end
	}

	@:arrayAccess
	public inline function get(i:Int):RealFloat {
		return this[i];
	}

	@:arrayAccess
	public inline function set(i:Int, v:RealFloat):RealFloat {
		this[i] = v;
		return v;
	}
}

@:forward(length)
private abstract FloatContainer(FloatContainerBase) {
	public inline function new(size:Int) {
		#if js
		this = new js.lib.Float64Array(size);
		#elseif sys
		this = new NativeContainer(size);
		#end
	}

	@:arrayAccess
	public inline function get(i:Int):Float {
		return this[i];
	}

	@:arrayAccess
	public inline function set(i:Int, v:Float):Float {
		this[i] = v;
		return v;
	}
}

private final class NativeContainerTools {
	public static inline function fill<T>(c:NativeContainer<T>, value:T):NativeContainer<T> {
		for (i in 0...c.length) {
			c[i] = value;
		}
		return c;
	}

	public static inline function subarray<T>(c:NativeContainer<T>, from:Int, to:Int):NativeContainer<T> {
		final container = new NativeContainer<T>(to - from);
		for (i in from...to) {
			container[i] = c[i];
		}
		return container;
	}
}
