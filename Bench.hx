package;

import delaunator.Delaunator;
import haxe.Timer;

final class Bench {
	static function main():Void
		new Bench();

	public function new() {
		#if js
		haxe.Log.trace = (a, ?b) -> untyped console.log(a);
		#elseif sys
		haxe.Log.trace = (a, ?b) -> Sys.println(a);
		#end

		final distributions = [uniform, gaussian, grid, degenerate];
		final names = ["uniform", "gaussian", "grid", "degenerate"];
		final counts = [20000, 100000, 200000, 500000, 1000000];

		for (i in 0...distributions.length) {
			final generate = distributions[i];
			trace(names[i]);

			// warmup
			Delaunator.from(generate(counts[0]));
			Delaunator.from(generate(counts[1]));

			for (c in counts) {
				final points = generate(c);
				final time = Timer.stamp();

				Delaunator.from(points);

				final ms = Timer.stamp() - time;
				trace('${fixedFloat(ms * 1000, 2)} ms');
			}
			trace("\n");
		}
	}

	function uniform(count:Int):Array<Array<Float>> {
		final points:Array<Array<Float>> = [];
		for (i in 0...count) {
			points.push([Math.random() * 1000, Math.random() * 1000]);
		}
		return points;
	}

	function grid(count:Int):Array<Array<Float>> {
		final points:Array<Array<Float>> = [];
		final size = Std.int(Math.sqrt(count));
		for (i in 0...size) {
			for (j in 0...size) {
				points.push([i, j]);
			}
		}
		return points;
	}

	function gaussian(count:Int):Array<Array<Float>> {
		final points:Array<Array<Float>> = [];
		for (i in 0...count) {
			points.push([pseudoNormal() * 1000, pseudoNormal() * 1000]);
		}
		return points;
	}

	function degenerate(count:Int):Array<Array<Float>> {
		final points:Array<Array<Float>> = [[0, 0]];
		for (i in 0...count) {
			final angle = 2 * Math.PI * i / count;
			points.push([10000000000 * Math.sin(angle), 10000000000 * Math.cos(angle)]);
		}
		return points;
	}

	function pseudoNormal():Float {
		final v = Math.random() + Math.random() + Math.random() + Math.random() + Math.random() + Math.random();
		return Math.min(0.5 * (v - 3) / 3, 1);
	}

	inline function fixedFloat(v:Float, precision:Int):Float {
		return Math.round(v * Math.pow(10, precision)) / Math.pow(10, precision);
	}
}
