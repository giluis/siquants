const std = @import("std");
const testing = std.testing;
const print = std.debug.print;
const utils = @import("utils.zig");
const SIOrders = @import("siorders.zig").SIOrders;
const quant = @import("quant.zig");
const QuantityKind = @import("quantkind.zig").QuantityKind; 
const Q = quant.Q;
const assert = std.debug.assert;
const Length = @import("quantkind.zig").Length;

fn Constant(comptime qk: QuantityKind, comptime unit: SIOrders, value: f64) Q(qk, .Unit){
    return Q(qk, unit ).init(value);
}

test "length * const" {
    const len = Length(.Milli).init(3.4);
    const planks = Constant(.Length, .Unit, 1.23e-10);
    const result = len.mult(planks) catch unreachable;

    assert(@TypeOf(result) == Q(.Area, .Nano));
    assert(equalUpToDigits(result.value, 0.000000000123,  9));

    print("\n {}", .{ @TypeOf(result)});
}

fn equalUpToDigits(a: f64, b: f64, digits: u32) bool {
    const factor = std.math.pow(f64, 10.0, @floatFromInt(digits));
    return @abs(a - b) * factor < 1.0;
}