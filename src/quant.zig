const std = @import("std");
const print = std.debug.print;
const fieldIndex = @import("utils.zig").fieldIndex;
const SIOrders = @import("siorders.zig").SIOrders;
const requireStructField = @import("utils.zig").require_struct_field;
const QuantityKind = @import("quantkind.zig").QuantityKind;
const Dimensions = @import("dimensions.zig").Dimensions;
const VALID_MULTS = @import("valid_mults.zig").VALID_MULTS;


pub fn Q(comptime qk: QuantityKind, comptime b: SIOrders, Self: type) type {
    return struct {
        const base: SIOrders = b;

        pub fn init(value: f64) Self {
            return Self {
                .value = value
            };
        }

        pub fn to_dimensions() Dimensions {
            return qk.to_dimensions();
        }

        pub fn as(self: Self, comptime result_unit: SIOrders) f64 {
            const mult_factor = comptime SIOrders.multiplication_factor(base, result_unit);
            return self.value * mult_factor;
        }

        pub fn from(comptime input_unit: SIOrders, value: f64) Self {
            const mult_factor = comptime SIOrders.multiplication_factor(input_unit, base);
            return Self{ .value = value * mult_factor };
        }


        pub fn get_base() SIOrders {
            return b;
        }

        pub fn get_qk() QuantityKind {
            return qk;
        }

        // TODO(luis.gil): add checked math for floats
        // TODO(luis.gil): verify that precision loss is acceptable after normalization
        // TODO(luis.gil): test this
        pub fn mult(self: Self, other: anytype) error{Overflow}! Q( QuantityKind.mult_quant_kinds_unchecked(Self.get_qk(), @TypeOf(other).get_qk()), b) {
            const new_qk_op = comptime QuantityKind.mult_quant_kinds(Self.get_qk(), @TypeOf(other).get_qk());
            const result_qk: QuantityKind = new_qk_op orelse @compileError("cannot mult " ++ @typeName(Self) ++ " with " ++ @typeName(@TypeOf(other)));

            const other_base = @TypeOf(other).get_base();
            const mult_factor = SIOrders.multiplication_factor(other_base, base);
            print("{}", .{b});

            return Q(result_qk,b) {
              .value = self.value * other.value * mult_factor,
            };
        }

        // TODO(luis.gil): handle overflows and the like
        pub fn scale(self: *Self, factor: f64) error{Overflow}!Self {
            return Q(qk, b).init(self.value * factor); 
        }

        // // TODO(luis.gil): add checked math for floats
        // // TODO(luis.gil): verify that precision loss is acceptable after normalization
        // // TODO(luis.gil): test this
        // pub fn div(self: Self, other: anytype) error{Overflow}!match_types_checked(Self, @TypeOf(other)) {
        //     const Ty = match_types(Self, @TypeOf(other)) orelse @compileError("cannot mult " ++ @typeName(Self) ++ " with " ++ @typeName(@TypeOf(other)));
        //     const other_base = @TypeOf(other).get_base();
        //     const multiplication_factor = SIOrders.multiplication_factor(other_base, base);

        //     return Ty{
        //         .value = self.value / (other.value * multiplication_factor),
        //     };
        // }

        pub fn to_value(self: Self) f64 {
            return self.value;
        }
    };
}


// /// Build a node type whose travel-depth is `degree`.
// pub fn Node(comptime degree: comptime_int) type {
//     // Leaf node: just a value and an `init`.
//     if (degree == 0) {
//         return struct {
//             value: u32,

//             pub fn init(v: u32) @This() {
//                 return .{ .value = v };
//             }
//         };
//     }

//     // Internal node: same payload + `travel`.
//     const Child = Node(degree - 1); // type for the next level down

//     return struct {
//         value: u32,

//         pub fn init(v: u32) @This() {
//             return .{ .value = v };
//         }

//         /// Move down one level, cloning the payload.
//         pub fn travel(self: @This()) Child {
//             return Child.init(self.value);
//         }
//     };
// }
