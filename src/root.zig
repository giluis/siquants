const std = @import("std");
const testing = std.testing;
const print = std.debug.print;
const utils = @import("utils.zig");
const removeNonUnique = utils.removeNonUnique;
const requireStructField = utils.require_struct_field;
const SIOrders = @import("siorders.zig").SIOrders;
const quant = @import("quant.zig");
const QuantityKind = quant.QuantityKind;
const Q = quant.Q;
const a = @import("valid_mults.zig").a;

test "length * time" {
    a();

    // const lenmili = Q(.Length, .Milli).init(2.0);
    // const areaunit = Q(.Area, .Unit).init(1.0);
    // const result_volume = lenmili.mult(areaunit) catch unreachable;
    // print("{}", result_volume.as(.Centi));
}
