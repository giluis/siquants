const std = @import("std");
const print = std.debug.print;
const fieldIndex = @import("utils.zig").fieldIndex;
const SIOrders = @import("siorders.zig").SIOrders;
const requireStructField = @import("utils.zig").require_struct_field;
const VALID_MULTS = @import("valid_mults.zig").VALID_MULTS;

pub const QuantityKind = enum {
    AbsorbedDoseRate,
    Action,
    AmountOfSubstance,
    AngularAcceleration,
    Area,
    AreaDensity,
    Capacitance,
    CatalyticActivityConcentration,
    ChemicalPotential,
    VolumeDensity,
    DoseEquivalent,
    ElectricCharge,
    ElectricChargeDensity,
    ElectricalConductance,
    ElectricalConductivity,
    ElectricCurrent,
    ElectricPotential,
    ElectricalResistance,
    ElectricalResistivity,
    Energy,
    EnergyDensity,
    Entropy,
    Frequency,
    HalfLife,
    Heat,
    HeatCapacity,
    HeatFluxDensity,
    Illuminance,
    Impedance,
    Inductance,
    Irradiance,
    Intensity,
    KineticEnergy,
    Length,
    LinearDensity,
    LuminousFlux,
    LuminousIntensity,
    MagneticFlux,
    Mass,
    MeanLifetime,
    MolarConcentration,
    MolarEnergy,
    MolarEntropy,
    MolarHeatCapacity,
    MomentOfInertia,
    OpticalPower,
    Permeability,
    Permittivity,
    PotentialEnergy,
    Power,
    Pressure,
    RadioActivity,
    RadiationDose,
    Radiance,
    RadiantIntensity,
    ReactionRate,
    Reluctance,
    SpecificEnergy,
    SpecificHeatCapacity,
    SpecificVolume,
    Speed,
    Spin,
    Stress,
    SurfaceTension,
    Temperature,
    ThermalConductance,
    ThermalConductivity,
    ThermalResistance,
    ThermalResistivity,
    Time,
    Viscosity,
    Volume,
    VolumetricFlowRate,
    Wavelength,
    Wavenumber,
    Work,
    YoungsModulus,
    SpringConstant,

    pub fn to_dimensions(self: QuantityKind) Dimensions {
        var d = Dimensions.zeros();
        return switch (self) {
            .AbsorbedDoseRate => d.setTime(-3).setLength(2),
            .Action => d.setTime(-1).setLength(2).setMass(1),
            .AmountOfSubstance => d.setAmountOfSubstance(1),
            .AngularAcceleration => d.setTime(-2),
            .Area => d.setLength(2),
            .AreaDensity => d.setLength(-2).setMass(1),
            .Capacitance => d.setTime(4).setLength(-2).setMass(-1).setElectricCurrent(2),
            .CatalyticActivityConcentration => d.setTime(-1).setLength(-3).setAmountOfSubstance(1),
            .ChemicalPotential => d.setTime(-2).setLength(2).setMass(1).setAmountOfSubstance(-1),
            .VolumeDensity => d.setLength(-3).setMass(1),
            .DoseEquivalent => d.setTime(-2).setLength(2),
            .ElectricCharge => d.setTime(1).setElectricCurrent(1),
            .ElectricChargeDensity => d.setTime(1).setLength(-3).setElectricCurrent(1),
            .ElectricalConductance => d.setTime(3).setLength(-2).setMass(-1).setElectricCurrent(2),
            .ElectricalConductivity => d.setTime(3).setLength(-3).setMass(-1).setElectricCurrent(2),
            .ElectricCurrent => d.setElectricCurrent(1),
            .ElectricPotential => d.setTime(-3).setLength(2).setMass(1).setElectricCurrent(-1),
            .ElectricalResistance => d.setTime(-3).setLength(2).setMass(1).setElectricCurrent(-2),
            .ElectricalResistivity => d.setTime(-3).setLength(3).setMass(1).setElectricCurrent(-2),
            .Energy => d.setTime(-2).setLength(2).setMass(1),
            .EnergyDensity => d.setTime(-2).setLength(-1).setMass(1),
            .Entropy => d.setTime(-2).setLength(2).setMass(1).setAbsoluteTemperature(-1),
            .Frequency => d.setTime(-1),
            .HalfLife => d.setTime(1),
            .Heat => d.setTime(-2).setLength(2).setMass(1),
            .HeatCapacity => d.setTime(-2).setLength(2).setMass(1).setAbsoluteTemperature(-1),
            .HeatFluxDensity => d.setTime(-3).setMass(1),
            .Illuminance => d.setLength(-2).setLuminousIntensity(1),
            .Impedance => d.setTime(-3).setLength(2).setMass(1).setElectricCurrent(-2),
            .Inductance => d.setTime(-2).setLength(2).setMass(1).setElectricCurrent(-2),
            .Irradiance => d.setTime(-3).setMass(1),
            .Intensity => d.setTime(-3).setMass(1),
            .KineticEnergy => d.setTime(-2).setLength(2).setMass(1),
            .Length => d.setLength(1),
            .LinearDensity => d.setLength(-1).setMass(1),
            .LuminousFlux => d.setLuminousIntensity(1),
            .LuminousIntensity => d.setLuminousIntensity(1),
            .MagneticFlux => d.setTime(-2).setLength(2).setMass(1).setElectricCurrent(-1),
            .Mass => d.setMass(1),
            .MeanLifetime => d.setTime(1),
            .MolarConcentration => d.setLength(-3).setAmountOfSubstance(1),
            .MolarEnergy => d.setTime(-2).setLength(2).setMass(1).setAmountOfSubstance(-1),
            .MolarEntropy => d.setTime(-2).setLength(2).setMass(1).setAbsoluteTemperature(-1).setAmountOfSubstance(-1),
            .MolarHeatCapacity => d.setTime(-2).setLength(2).setMass(1).setAbsoluteTemperature(-1).setAmountOfSubstance(-1),
            .MomentOfInertia => d.setLength(2).setMass(1),
            .OpticalPower => d.setLength(-1),
            .Permeability => d.setTime(-2).setLength(1).setMass(1).setElectricCurrent(-2),
            .Permittivity => d.setTime(4).setLength(-3).setMass(-1).setElectricCurrent(2),
            .PotentialEnergy => d.setTime(-2).setLength(2).setMass(1),
            .Power => d.setTime(-3).setLength(2).setMass(1),
            .Pressure => d.setTime(-2).setLength(-1).setMass(1),
            .RadioActivity => d.setTime(-1),
            .RadiationDose => d.setTime(-2).setLength(2),
            .Radiance => d.setTime(-3).setMass(1),
            .RadiantIntensity => d.setTime(-3).setLength(2).setMass(1),
            .ReactionRate => d.setTime(-1).setLength(-3).setAmountOfSubstance(1),
            .Reluctance => d.setTime(2).setLength(-2).setMass(-1).setElectricCurrent(2),
            .SpecificEnergy => d.setTime(-2).setLength(2),
            .SpecificHeatCapacity => d.setTime(-2).setLength(2).setAbsoluteTemperature(-1),
            .SpecificVolume => d.setLength(3).setMass(-1),
            .Speed => d.setTime(-1).setLength(1),
            .Spin => d.setTime(-1).setLength(2).setMass(1),
            .Stress => d.setTime(-2).setLength(-1).setMass(1),
            .SurfaceTension => d.setTime(-2).setMass(1),
            .Time => d.setTime(1),
            .Temperature => d.setAbsoluteTemperature(1),
            .ThermalConductance => d.setTime(-3).setLength(2).setMass(1).setAbsoluteTemperature(-1),
            .ThermalConductivity => d.setTime(-3).setLength(1).setMass(1).setAbsoluteTemperature(-1),
            .ThermalResistance => d.setTime(3).setLength(-2).setMass(-1).setAbsoluteTemperature(1),
            .ThermalResistivity => d.setTime(3).setLength(-1).setMass(-1).setAbsoluteTemperature(1),
            .Viscosity => d.setTime(-1).setLength(-1).setMass(1),
            .Volume => d.setLength(3),
            .VolumetricFlowRate => d.setTime(-1).setLength(3),
            .Wavelength => d.setLength(1),
            .Wavenumber => d.setLength(-1),
            .Work => d.setTime(-2).setLength(2).setMass(1),
            .YoungsModulus => d.setTime(-2).setLength(-1).setMass(1),
            .SpringConstant => d.setTime(-2).setMass(1),
        };
    }

    pub fn mult_quant_kinds(first_qk: @This(), second_qk: @This()) ?@This() {
        const slice= &VALID_MULTS[@intFromEnum(first_qk)];
        for (&slice) |item| {
            if (item[0] == second_qk) {
                return item[1];
            }
        }
        return null;
    }
};

pub fn Q(comptime qk: QuantityKind, comptime b: SIOrders) type {
    // requireStructField(Self, "value", f64);

    return struct {
        const Self = @This();
        const base: SIOrders = b;

        value: f64,

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
        pub fn mult(self: Self, other: anytype) error{Overflow}! Self {
            const new_qk_op = comptime QuantityKind.mult_quant_kinds(Self.get_qk(), @TypeOf(other).get_qk());
            const result_qk = undefined;
            if (new_qk_op) |new_qk|{
                result_qk = new_qk;
            } else {
                @compileError("cannot mult " ++ @typeName(Self) ++ " with " ++ @typeName(@TypeOf(other)));
            }

            const other_base = @TypeOf(other).get_base();
            const mult_factor = SIOrders.multiplication_factor(other_base, base);

            return Q(new_qk_op,b) {
                .value = self.value * other.value * mult_factor,
            };
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

/// Each 4 bit represents an exponent from -7 to 8
/// Each physical measure can be described by combinations of these exponentes
/// E.g acceleration = d / t^2 i.e so it would be represented as
/// Dimensions {
///     .time: -2
///     .length: 1
///     ... all other fields to zero
/// }
//
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
        print("\n=={s}", .{input});

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
