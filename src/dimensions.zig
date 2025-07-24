const std = @import("std");
const fieldIndex = @import("utils.zig").fieldIndex;
const QuantityKind = @import("quantkind.zig").QuantityKind;

/// This type represents the dimensional algebraic representation of each quantity kind as 
/// 
/// Each 4 bit represents an exponent from -7 to 8
/// Each physical measure can be described by combinations of these exponentes
/// E.g acceleration = d / t^2 i.e so it would be represented as
/// Dimensions {
///     .time: -2
///     .length: 1
///     ... all other fields to zero
/// }
/// 
/// SAFETY: TODO add docs for this
pub const Dimensions = packed struct {
    mass: i4,
    length: i4,
    time: i4,
    electric_current: i4,
    absolute_temperature: i4,
    amount_of_substance: i4,
    luminous_intensity: i4,
    _padding: i4,

    const Self = @This();

    pub fn addChecked(comptime T: type, a: T, b: T) error{Overflow}!T {
        const result = @addWithOverflow(a, b);
        if (result.@"1" == 1) {
            return error.Overflow;
        }
        return result.@"0";
    }

    pub fn add_dimensions(self: Self, other: Self) error{Overflow}!Self {
        return Self{
            .time = try Self.addChecked(i4, self.time, other.time),
            .length = try Self.addChecked(i4, self.length, other.length),
            .electric_current = try Self.addChecked(i4, self.electric_current, other.electric_current),
            .absolute_temperature = try Self.addChecked(i4, self.absolute_temperature, other.absolute_temperature),
            .luminous_intensity = try Self.addChecked(i4, self.luminous_intensity, other.luminous_intensity),
            .amount_of_substance = try Self.addChecked(i4, self.amount_of_substance, other.amount_of_substance),
            .mass = try Self.addChecked(i4, self.mass, other.mass),
            ._padding = 0,
        };
    }

    pub fn iter_fields(self: Self) [@typeInfo(Self).@"struct".fields.len - 1]i4 {
        const field_count = @typeInfo(Self).@"struct".fields.len - 1;
        var result: [field_count]i4 = undefined;

        const fields = @typeInfo(Self).@"struct".fields;

        comptime var i: usize = 0;
        inline while (i < field_count) : (i += 1) {
            const name = fields[i].name;
            result[i] = @field(self, name);
        }

        return result;
    }

    pub fn setTime(self: Self, time_exponent: i4) Self {
        var other = self;
        other.time = time_exponent;
        return other;
    }

    pub fn setLength(self: Self, length_exponent: i4) Self {
        var other = self;
        other.length = length_exponent;
        return other;
    }

    pub fn setMass(self: Self, mass_exponent: i4) Self {
        var other = self;
        other.mass = mass_exponent;
        return other;
    }

    pub fn setAmountOfSubstance(self: Self, aos_exponent: i4) Self {
        var other = self;
        other.amount_of_substance = aos_exponent;
        return other;
    }

    pub fn setAbsoluteTemperature(self: Self, absolute_temp: i4) Self {
        var other = self;
        other.absolute_temperature = absolute_temp;
        return other;
    }

    pub fn setElectricCurrent(self: Self, electric_current: i4) Self {
        var other = self;
        other.electric_current = electric_current;
        return other;
    }

    pub fn setLuminousIntensity(self: Self, luminous_intensity: i4) Self {
        var other = self;
        other.luminous_intensity = luminous_intensity;
        return other;
    }

    const FromStrErr = error{DimensionLess};

    fn write(input: []u8, at_arg: usize, char: u8, value: i4) usize {
        var at = at_arg;
        if (value == 0) {
            return at;
        }
        if (value == 1) {
            input[at] = char;
            return at + 1;
        }

        input[at] = char;
        at += 1;
        if (value < 0) {
            input[at] = '-';
            at += 1;
        }
        input[at] = @as(u8, @abs(value)) + '0';
        return at + 1;
    }

    pub fn to_string_validation(self: *const Self, allocator: std.mem.Allocator) []u8 {
        const result = allocator.alloc(u8, 80) catch unreachable;
        var at: usize = 0;
        inline for (@typeInfo(Self).@"struct".fields) |field| {
            at = write(result, at, std.ascii.toUpper(field.name[0]), @field(self.*, field.name));
        }

        return result[0..at];
    }

    pub fn from_str(input: []const u8) FromStrErr!Self {
        var dims: [7]i4 = .{0} ** 7;

        const length_idx = comptime fieldIndex(Dimensions, "length");
        const time_idx = comptime fieldIndex(Dimensions, "time");
        const mass_idx = comptime fieldIndex(Dimensions, "mass");
        const electric_current_idx = comptime fieldIndex(Dimensions, "electric_current");
        const absolute_temperature_idx = comptime fieldIndex(Dimensions, "absolute_temperature");
        const amount_of_substance_idx = comptime fieldIndex(Dimensions, "amount_of_substance");
        const luminous_intensity_idx = comptime fieldIndex(Dimensions, "luminous_intensity");

        var result = std.mem.tokenizeScalar(u8, input, ' ');

        while (result.next()) |dim_str| {
            const dim_idx = switch (dim_str[0]) {
                'L' => length_idx,
                'T' => time_idx,
                'I' => electric_current_idx,
                'O' => absolute_temperature_idx,
                'M' => mass_idx,
                'N' => amount_of_substance_idx,
                'J' => luminous_intensity_idx,
                '1' => return FromStrErr.DimensionLess,
                else => |v| {
                    std.debug.print("v: {c}", .{v});
                    unreachable;
                },
            };

            var dim: i8 = undefined;

            if (dim_str.len == 1) {
                dim = 1;
            } else {
                const sign: i8 = if (dim_str.len == 2) 1 else -1;
                const exp_idx: usize = if (sign == -1) 2 else 1;
                const tempu8: i8 = @intCast(dim_str[exp_idx] - '0');
                dim = tempu8 * sign;
            }
            dims[dim_idx] = @intCast(dim);
        }

        return Self{
            .time = dims[0],
            .length = dims[1],
            .mass = dims[2],
            .electric_current = dims[3],
            .absolute_temperature = dims[4],
            .amount_of_substance = dims[5],
            .luminous_intensity = dims[6],
            ._padding = 0,
        };
    }

    pub fn zeros() Self {
        return Dimensions{
            .time = 0,
            .length = 0,
            .mass = 0,
            .electric_current = 0,
            .absolute_temperature = 0,
            .amount_of_substance = 0,
            .luminous_intensity = 0,
            ._padding = 0,
        };
    }

    fn toUnit(self: Self) QuantityKind {
        switch (self) {
            .{
                .time = 0,

                .length = 1,
                .mass = 0,
                .electric_current = 0,
                .absolute_temperature = 0,
                .amount_of_substance = 0,
                .luminous_intensity = 0,
                ._padding = 0,
            } => return QuantityKind.Length,
            else => return QuantityKind.AbsorbedDose,
        }
    }
};