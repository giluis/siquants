const Dimensions = @import("dimensions.zig").Dimensions;
const VALID_MULTS = @import("valid_mults.zig").VALID_MULTS;
const SIOrders = @import("siorders.zig").SIOrders;
const Q = @import("quant.zig").Q;

pub const QuantityKind = enum {
    AbsorbedDoseRate,
    Acceleration,
    Action,
    AmountOfSubstance,
    AngularAcceleration,
    AngularMomemtum,
    Area,
    AreaDensity,
    Capacitance,
    CatalyticActivityConcentration,
    ChemicalPotential,
    Compressibility,
    VolumeDensity,
    DoseEquivalent,
    Density,
    DynamicViscosity,
    ElectricCharge,
    ElectricChargeDensity,
    ElectricField,
    ElectricalConductance,
    ElectricalConductivity,
    ElectricCurrent,
    ElectricPotential,
    ElectricalResistance,
    ElectricalResistivity,
    Energy,
    EnergyDensity,
    EnergyPerMass,
    Entropy,
    Frequency,
    Force,
    HalfLife,
    // todo: should I include heat?
    Heat,
    HeatCapacity,
    HeatFluxDensity,
    Illuminance,
    Impedance,
    Impulse,
    Inductance,
    Irradiance,
    Intensity,
    // todo: should I add kinectic energy?
    KineticEnergy,
    LatentHeat,
    Length,
    LinearDensity,
    LuminousFlux,
    LuminousIntensity,
    MagneticFlux,
    Mass,
    MagneticField,
    // todo: should I keep mean lifetime?
    // MeanLifetime,
    MolarConcentration,
    MolarEnergy,
    MolarEntropy,
    MolarHeatCapacity,
    MomentOfInertia,
    Momemtum,
    NumberOfCycles,
    OpticalPower,
    Permeability,
    Permittivity,
    // todo: should I add potential energy?
    // PotentialEnergy,
    Power,
    Pressure,
    RadioActivity,
    RadiationDose,
    Radiance,
    RadiantIntensity,
    ReactionRate,
    Reluctance,
    SolidAngle,
    SpecificEnergy,
    SpecificHeatCapacity,
    SpecificVolume,
    Speed,
    Spin,
    Stress,
    Strain,
    SurfaceTension,
    Temperature,
    TemperatureDiff,
    ThermalConductance,
    ThermalConductivity,
    ThermalResistance,
    ThermalResistivity,
    Time,
    Velocity,
    Viscosity,
    Volume,
    VolumetricFlowRate,
    Wavelength,
    Wavenumber,
    // todo: should I have work?
    Work,
    YoungsModulous,
    SpringConstant,

    UNKNOWN,

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

    pub fn mult_quant_kinds(comptime first_qk: @This(), comptime second_qk: @This()) ?@This() {
        const a = @tagName(first_qk);
        const slice = @field(VALID_MULTS, a);
        for (&slice) |item| {
            if (item[0] == second_qk) {
                return item[1];
            }
        }
        return null;
    }

    pub fn mult_quant_kinds_unchecked(comptime first_qk: @This(), comptime second_qk: @This()) @This() {
        return mult_quant_kinds(first_qk, second_qk) orelse .UNKNOWN;
    }
};

pub fn AbsorbedDoseRate(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.AbsorbedDoseRate, base, Self);
    };
}

pub fn Acceleration(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Acceleration, base, Self);
    };
}

pub fn Action(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Action, base, Self);
    };
}

pub fn AmountOfSubstance(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.AmountOfSubstance, base, Self);
    };
}

pub fn AngularAcceleration(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.AngularAcceleration, base, Self);
    };
}

pub fn AngularMomemtum(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.AngularMomemtum, base, Self);
    };
}

pub fn Area(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Area, base, Self);
    };
}

pub fn AreaDensity(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.AreaDensity, base, Self);
    };
}

pub fn Capacitance(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Capacitance, base, Self);
    };
}

pub fn CatalyticActivityConcentration(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.CatalyticActivityConcentration, base, Self);
    };
}

pub fn ChemicalPotential(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.ChemicalPotential, base, Self);
    };
}

pub fn Compressibility(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Compressibility, base, Self);
    };
}

pub fn VolumeDensity(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.VolumeDensity, base, Self);
    };
}

pub fn DoseEquivalent(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.DoseEquivalent, base, Self);
    };
}

pub fn Density(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Density, base, Self);
    };
}

pub fn DynamicViscosity(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.DynamicViscosity, base, Self);
    };
}

pub fn ElectricCharge(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.ElectricCharge, base, Self);
    };
}

pub fn ElectricChargeDensity(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.ElectricChargeDensity, base, Self);
    };
}

pub fn ElectricField(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.ElectricField, base, Self);
    };
}

pub fn ElectricalConductance(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.ElectricalConductance, base, Self);
    };
}

pub fn ElectricalConductivity(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.ElectricalConductivity, base, Self);
    };
}

pub fn ElectricCurrent(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.ElectricCurrent, base, Self);
    };
}

pub fn ElectricPotential(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.ElectricPotential, base, Self);
    };
}

pub fn ElectricalResistance(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.ElectricalResistance, base, Self);
    };
}

pub fn ElectricalResistivity(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.ElectricalResistivity, base, Self);
    };
}

pub fn Energy(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Energy, base, Self);
    };
}

pub fn EnergyDensity(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.EnergyDensity, base, Self);
    };
}

pub fn EnergyPerMass(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.EnergyPerMass, base, Self);
    };
}

pub fn Entropy(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Entropy, base, Self);
    };
}

pub fn Frequency(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Frequency, base, Self);
    };
}

pub fn Force(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Force, base, Self);
    };
}

pub fn HalfLife(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.HalfLife, base, Self);
    };
}

pub fn Heat(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Heat, base, Self);
    };
}

pub fn HeatCapacity(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.HeatCapacity, base, Self);
    };
}

pub fn HeatFluxDensity(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.HeatFluxDensity, base, Self);
    };
}

pub fn Illuminance(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Illuminance, base, Self);
    };
}

pub fn Impedance(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Impedance, base, Self);
    };
}

pub fn Impulse(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Impulse, base, Self);
    };
}

pub fn Inductance(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Inductance, base, Self);
    };
}

pub fn Irradiance(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Irradiance, base, Self);
    };
}

pub fn Intensity(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Intensity, base, Self);
    };
}

pub fn KineticEnergy(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.KineticEnergy, base, Self);
    };
}

pub fn LatentHeat(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.LatentHeat, base, Self);
    };
}

pub fn Length(comptime base: SIOrders) type {
    return struct {
        const Self = @This();         // alias the enclosing type
        value: f64,                   // storage

        // pull in all of Q’s methods (init, as, from, mult, etc.)
        usingnamespace Q(.Length, base, Self);
    };
}


pub fn LinearDensity(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.LinearDensity, base, Self);
    };
}

pub fn LuminousFlux(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.LuminousFlux, base, Self);
    };
}

pub fn LuminousIntensity(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.LuminousIntensity, base, Self);
    };
}

pub fn MagneticFlux(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.MagneticFlux, base, Self);
    };
}

pub fn Mass(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Mass, base, Self);
    };
}

pub fn MagneticField(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.MagneticField, base, Self);
    };
}

pub fn MolarConcentration(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.MolarConcentration, base, Self);
    };
}

pub fn MolarEnergy(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.MolarEnergy, base, Self);
    };
}

pub fn MolarEntropy(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.MolarEntropy, base, Self);
    };
}

pub fn MolarHeatCapacity(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.MolarHeatCapacity, base, Self);
    };
}

pub fn MomentOfInertia(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.MomentOfInertia, base, Self);
    };
}

pub fn Momemtum(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Momemtum, base, Self);
    };
}

pub fn NumberOfCycles(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.NumberOfCycles, base, Self);
    };
}

pub fn OpticalPower(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.OpticalPower, base, Self);
    };
}

pub fn Permeability(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Permeability, base, Self);
    };
}

pub fn Permittivity(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Permittivity, base, Self);
    };
}

pub fn Power(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Power, base, Self);
    };
}

pub fn Pressure(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Pressure, base, Self);
    };
}

pub fn RadioActivity(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.RadioActivity, base, Self);
    };
}

pub fn RadiationDose(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.RadiationDose, base, Self);
    };
}

pub fn Radiance(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Radiance, base, Self);
    };
}

pub fn RadiantIntensity(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.RadiantIntensity, base, Self);
    };
}

pub fn ReactionRate(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.ReactionRate, base, Self);
    };
}

pub fn Reluctance(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Reluctance, base, Self);
    };
}

pub fn SolidAngle(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.SolidAngle, base, Self);
    };
}

pub fn SpecificEnergy(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.SpecificEnergy, base, Self);
    };
}

pub fn SpecificHeatCapacity(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.SpecificHeatCapacity, base, Self);
    };
}

pub fn SpecificVolume(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.SpecificVolume, base, Self);
    };
}

pub fn Speed(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Speed, base, Self);
    };
}

pub fn Spin(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Spin, base, Self);
    };
}

pub fn Stress(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Stress, base, Self);
    };
}

pub fn Strain(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Strain, base, Self);
    };
}

pub fn SurfaceTension(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.SurfaceTension, base, Self);
    };
}

pub fn Temperature(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Temperature, base, Self);
    };
}

pub fn TemperatureDiff(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.TemperatureDiff, base, Self);
    };
}

pub fn ThermalConductance(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.ThermalConductance, base, Self);
    };
}

pub fn ThermalConductivity(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.ThermalConductivity, base, Self);
    };
}

pub fn ThermalResistance(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.ThermalResistance, base, Self);
    };
}

pub fn ThermalResistivity(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.ThermalResistivity, base, Self);
    };
}

pub fn Time(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Time, base, Self);
    };
}

pub fn Velocity(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Velocity, base, Self);
    };
}

pub fn Viscosity(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Viscosity, base, Self);
    };
}

pub fn Volume(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Volume, base, Self);
    };
}

pub fn VolumetricFlowRate(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.VolumetricFlowRate, base, Self);
    };
}

pub fn Wavelength(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Wavelength, base, Self);
    };
}

pub fn Wavenumber(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Wavenumber, base, Self);
    };
}

pub fn Work(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.Work, base, Self);
    };
}

pub fn YoungsModulous(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.YoungsModulous, base, Self);
    };
}

pub fn SpringConstant(comptime base: SIOrders) type {
    return struct {
        const Self = @This();
        value: f64,
        usingnamespace Q(.SpringConstant, base, Self);
    };
}
