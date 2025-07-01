const POWERS_OF_10: [21]f64 = .{
    1e-24,
    1e-21,
    1e-18,
    1e-15,
    1e-12,
    1e-9,
    1e-6,
    1e-3,
    1e-2,
    1e-1,
    1e0,
    1e1,
    1e2,
    1e3,
    1e6,
    1e9,
    1e12,
    1e15,
    1e18,
    1e21,
    1e24,
};

pub const SIOrders = enum {
    Yocto,
    Zepto,
    Atto,
    Femto,
    Pico,
    Nano,
    Micro,
    Milli,
    Centi,
    Deci,
    Unit,
    Deca,
    Hecta,
    Kilo,
    Mega,
    Giga,
    Tera,
    Peta,
    Exa,
    Zetta,
    Yotta,

    // TODO(luis.gil): test this
    pub fn multiplication_factor(from: @This(), to: @This()) f64 {
        const int_from: i8 = @intFromEnum(from);
        const int_to: i8 = @intFromEnum(to);
        const diff: i8 = int_to - int_from;

        const enum_size = @typeInfo(@This()).@"enum".fields.len;
        const make_it_pos = @as(i8, enum_size - @intFromEnum(@This().Unit));
        const pos_diff: u8 = @intCast(diff + make_it_pos);

        const idx: usize = @intCast(pos_diff);
        return POWERS_OF_10[idx];
    }

    pub fn asUnit(from: f64) f64 {
        const mult_factor = SIOrders.multiplication_factor(from, .Unit);
        return from * mult_factor;
    }
};
