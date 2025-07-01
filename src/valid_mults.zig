const QuantityKind = @import("quant.zig").QuantityKind;
const std = @import("std");
const MultPair: type = [2]QuantityKind;

pub fn a() void {
    const mult_fields =@typeInfo(VALID_MULTS).@"struct".fields;
    const qk_fields =@typeInfo(QuantityKind).@"struct".fields;
    std.debug.assert(mult_fields.len == qk_fields.len);
    @compileLog("...");

    for(mult_fields, qk_fields) |m,q| {
        const lower_m = [0] ** 100;
        const lower_q = [0] ** 100;
        std.ascii.lowerString(lower_m, m.name);
        std.ascii.lowerString(lower_q, q.name);

        std.debug.assert(std.mem.eql(u8, lower_m, lower_q)) ;
    }
}

// TODO: check this is correct
pub const VALID_MULTS= struct {
    absorbedDoseRate: [6]MultPair= .{ .{ .Mass, .Power }, .{ .HalfLife, .RadiationDose } },
    action: [6]MultPair= .{.{ .Frequency, .Energy }},
    amountOfSubstance: [6]MultPair= .{
        .{ .ChemicalPotential, .Energy },
        .{ .MolarEnergy, .Energy },
        .{ .MolarEntropy, .Entropy },
        .{ .MolarHeatCapacity, .HeatCapacity },
    },
    angularAcceleration: [6]MultPair= .{},
    accelaration: [6]MultPair= .{
        .{ .Mass, .Force },
        .{ .Time, .Velocity },
        // Add Energy per mass
        .{ .Length, .EnergyPerMass },
    },
    area: [6]MultPair= .{
        .{ .Length, .Volume },
        .{ .Pressure, .Force },
        .{ .Illuminance, .LuminousFlux },
        .{ .MagneticField, .MagneticFlux },
        .{ .SurfaceTension, .Work },
    },
    areaDensity: [6]MultPair= .{.{ .Area, .Mass }},
    capacitance: [6]MultPair= .{
        .{ .ElectricPotential, .EletricCharge },
        // todo: consider Voltage^2 here. Also, does voltage^2 exist here?
        .{.ElectricPotentialSquared, .EletricCharge},
    },
    // todo: add Vector
    magneticField: [6]MultPair= .{
        .{ .Area, .MagneticFlux },
    },
    catalyticActivityConcentration: [6]MultPair= .{},
    chemicalPotential: [6]MultPair= .{
        .{ .MolarConcentration, .EnergyDensity },
    },
    volumeDensity: [6]MultPair= .{
        .{ .Volume, .Mass },
    },
    doseEquivalent: [6]MultPair= .{.{ .Mass, .Energy }},
    electricCharge: [6]MultPair= .{
        // todo: consider voltage herecross product of velocity and idk what else, need help here
        .{ .ElectricPotential, .Energy },
        // todo: add electric field
        .{ .ElectricField, .Force },
        // todo: this is the cross product of velocity and idk what else, need help herecross product of velocity and idk what else, need help here
        .{ .MagneticField, .Force },
    },
    electricChargeDensity: [6]MultPair= .{},
    electricalConductance: [6]MultPair= .{
        .{ .ElectricPotential, .ElectricCurrent },
    },
    electricalConductivity: [6]MultPair= .{},
    electricCurrent: [6]MultPair= .{
        .{ .Time, .ElectricCharge },
        // todo: consider voltage here
        .{ .ElectricResistance, .ElectricPotential },
        // todo: consider voltage here
        .{ .ElectricPotential, .Power },
    },
    electricField: [6]MultPair= .{
        .{ .Charge, .Force },
        // todo: consider electric potential
        .{ .Length, .Voltage },
    },
    // todo: consider voltage here
    electricPotential: [6]MultPair= .{
        .{ .ElectricCurrent, .Power },
        .{ .ElectricCharge, .Energy },
        .{ .ElectricalConductance, .ElectricCurrent },
        .{ .Time, .MagneticFlux },
    },
    electricalResistance: [6]MultPair= .{
        .{ .ElectricCurrent, .ElectricPotential },
        .{ .Length, .ElectricalResistivity },
        // todo: add a name for this
        .{.ElectricCurrentSquared, .ElectricPotential},
    },
    electricalResistivity: [6]MultPair= .{},
    energy: [6]MultPair= .{ .{ .Frequency, .Power }, .{ .Time, .Action } },
    energyDensity: [6]MultPair= .{.{ .Volume, .Energy }},
    entropy: [6]MultPair= .{
    // todo: should I have heat
    .{ .Temperature, .Heat }},
    frequency: [6]MultPair= .{
        .{ .PlancksConstant, .Energy },
        // todo: add number of cycles
        .{ .Time, .NumberOfCycles },
    },
    // todo: Vector
    force: [6]MultPair= .{
        .{ .Time, .Impulse },
        // todo: Should I separate length and distance
        // Applicable for Torque as well
        .{ .Length, .Energy },
        .{ .Velocity, .Power },
    },
    halfLife: [6]MultPair= .{.{ .AbsorbedDoseRate, .RadiationDose }},
    // heat = ,: [6]MultPair=
    heatCapacity: [6]MultPair= .{
    // todo: add heat?
        .{ .Temperature, .Heat }},
    heatFluxDensity: [6]MultPair= .{},
    illuminance: [6]MultPair= .{.{ .Area, .LuminousFlux }},
    impedance: [6]MultPair= .{.{ .ElectricCurrent, .ElectricPotential }},
    inductance: [6]MultPair= .{
        .{ .ElectricCurrent, .MagneticFlux },
        // todo: not sure if should add current squared
        .{.ElectricCurrentSquared, .MagneticFlux}
    },
    irradiance: [6]MultPair= .{.{ .Area, .Power }},
    intensity: [6]MultPair= .{.{ .Area, .Power }},
    // kineticEnergy = ,: [6]MultPair=
    length: [6]MultPair= .{ .{ .Length, .Area }, .{ .Area, .Volume }, .{ .Force, .Energy }, .{ .Momentum, .AngularMomemtum } },
    linearDensity: [6]MultPair= .{
        .{ .Length, .Mass },
    },
    luminousFlux: [6]MultPair= .{.{ .Time, .Energy }},
    luminousIntensity: [6]MultPair= .{
        // todo: add solid angle
        .{ .SolidAngle, .LuminousFlux },
    },
    magneticFlux: [6]MultPair= .{
        // todo: Inverse here
        // .{.InverseInductance, .Current}
    },
    mass: [6]MultPair= .{
        .{ .Acceleration, .Force },
        .{ .Velocity, .Momemtum },
        .{ .LatentHeat, .Energy },
    },
    // meanLifetime = ,: [6]MultPair=
    momemtum: [6]MultPair= .{
        // Add angular momentum
        .{ .Length, .AngularMomemtum },
        .{ .Velocity, .Energy },
    },
    molarConcentration: [6]MultPair= .{
        .{ .Volume, .AmountOfSubstance },
        .{ .ChemicalPotential, .EnergyDensity },
        .{ .MolarEnergy, .EnergyDensity },
    },
    molarEnergy: [6]MultPair= .{
        .{ .AmountOfSubstance, .Energy },
        .{ .MolarConcentration, .EnergyDensity },
    },
    molarEntropy: [6]MultPair= .{
        .{ .AmountOfSubstance, .Entropy },
    },
    molarHeatCapacity: [6]MultPair= .{
        .{ .AmountOfSubstance, .HeatCapacity },
    },
    momentOfInertia: [6]MultPair= .{},
    opticalPower: [6]MultPair= .{},
    permeability: [6]MultPair= .{},
    permittivity: [6]MultPair= .{},
    // TODO(luis.gil): add this
    // potentialEnergy: [6]MultPair= .{},
    power: [6]MultPair= .{
        .{ .Frequency, .Energy },
    },
    // todo: should i descriminate between energy, work and heat?
    pressure: [6]MultPair= .{
        .{ .Volume, .Energy },
        .{ .Area, .Force },
        // todo: check that dynamic viscosity exists
        .{ .Time, .DynamicViscosity },
        .{ .Length, .SurfaceTension },
        // todo: strain, compressibility
        .{ .Compressibility, .Strain },
    },
    radioActivity: [6]MultPair= .{},
    radiationDose: [6]MultPair= .{.{ .Mass, .Energy }},
    radiance: [6]MultPair= .{
        .{ .Area, .RadiantIntensity },
        .{ .SolidAngle, .Irradiance },
    },
    radiantIntensity: [6]MultPair= .{.{ .SolidAngle, .Power }},
    reactionRate: [6]MultPair= .{},
    reluctance: [6]MultPair= .{
        .{ .MagneticFlux, .ElectricCurrent },
    },
    specificEnergy: [6]MultPair= .{
        .{ .Mass, .Energy },
    },
    specificHeatCapacity: [6]MultPair= .{
        .{ .Mass, .HeatCapacity },
    },
    specificVolume: [6]MultPair= .{ .Mass, .Volume },
    speed: [6]MultPair= .{},
    spin: [6]MultPair= .{.{ .Action, .AngularMomentum }},
    stress: [6]MultPair= .{
        .{ .Volume, .Energy },
        .{ .Area, .Force },
        .{ .Strain, .EnergyDensity },
        .{ .Time, .Viscosity },
    },
    // todo: check strain exists and what is it
    strain: [6]MultPair= .{
        .{ .YoungsModulous, .Stress },
    },
    surfaceTension: [6]MultPair= .{
        .{ .Length, .Force },
        .{ .Area, .Energy },
    },
    temperature: [6]MultPair= .{},
    thermalConductance: [6]MultPair= .{
    // todo: Add temperature diff
    .{ .TemperatureDiff, .Power }},
    thermalConductivity: [6]MultPair= .{
        .{ .Length, .ThermalConductance },
    },
    thermalResistance: [6]MultPair= .{ .{ .Power, .TemperatureDiff }, .{ .Length, .ThermalResistivity } },
    thermalResistivity: [6]MultPair= .{},
    time: [6]MultPair= .{
        // Add impulse
        .{ .Force, .Impulse },
        .{ .Power, .Energy },
        .{ .Current, .EletricCharge },
        .{ .Acceleration, .Velocity },
        .{ .Frequency, .NumberOfCycles },
    },
    viscosity: [6]MultPair= .{
        // TODO: this is a division actually
        // .{.InverseTime, .Stress},
    },
    velocity: [6]MultPair= .{
        // Check displacement
        .{ .Time, .Length },
        .{ .Mass, .Momentum },
        .{ .Force, .Power },
    },
    volume: [6]MultPair= .{
        .{ .Pressure, .Word },
        .{ .Density, .Mass },
    },
    volumetricFlowRate: [6]MultPair= .{},
    wavelength: [6]MultPair= .{
        .{ .Frequency, .Speed },
    },
    wavenumber: [6]MultPair= .{},
    youngsModulus: [6]MultPair= .{
        .{ .Strain, .Stress },
    },
    springConstant: [6]MultPair= .{.{ .Length, .Force }},
    // todo: add solid angle
    solidAngle: [6]MultPair= .{},
};

    // // .AbsorbedDoseRate
    // // .Action
    // .{.{ .Frequency, .Energy }},
    // // .AmountOfSubstance
    // .{
    //     .{ .ChemicalPotential, .Energy },
    //     .{ .MolarEnergy, .Energy },
    //     .{ .MolarEntropy, .Entropy },
    //     .{ .MolarHeatCapacity, .HeatCapacity },
    // },
    // // .AngularAcceleration
    // // .Accelaration
    // // .Area
    // // .AreaDensity
    // // .Capacitance
    // // .MagneticField
    // // .CatalyticActivityConcentration
    // // .ChemicalPotential
    // // .VolumeDensity
    // // .DoseEquivalent
    // // .ElectricCharge
    // // .ElectricChargeDensity
    // // .ElectricalConductance
    // // .ElectricalConductivity
    // ,
    // // .ElectricCurrent
    // ,
    // // .ElectricField
    // ,
    // // todo: consider voltage here
    // // .ElectricPotential
    // ,
    // // .ElectricalResistance
    // ,
    // // .ElectricalResistivity
    // ,
    // // .Energy
    // ,
    // // .EnergyDensity
    // ,
    // // .Entropy
    // // todo: should I have heat
    // // .Frequency
    // ,
    // // todo: Vector
    // // .Force
    // .{
    //     .{ .Time, .Impulse },
    //     // todo: Should I separate length and distance
    //     // Applicable for Torque as well
    //     .{ .Length, .Energy },
    //     .{ .Velocity, .Power },
    // },
    // // .HalfLife
    // .{.{ .AbsorbedDoseRate, .RadiationDose }},
    // // .Heat = ,
    // // .HeatCapacity
    // .{
    // // todo: add heat?
    // .{ .Temperature, .Heat }},
    // // .HeatFluxDensity
    // .{},
    // // .Illuminance
    // .{.{ .Area, .LuminousFlux }},
    // // .Impedance
    // .{.{ .ElectricCurrent, .ElectricPotential }},
    // // .Inductance
    // .{
    //     .{ .ElectricCurrent, .MagneticFlux },
    //     // todo: not sure if should add current squared
    //     // .{.ElectricCurrentSquared, .MagneticFlux},
    // },
    // // .Irradiance
    // .{.{ .Area, .Power }},
    // // .Intensity
    // .{.{ .Area, .Power }},
    // // .KineticEnergy = ,
    // // .Length
    // .{ .{ .Length, .Area }, .{ .Area, .Volume }, .{ .Force, .Energy }, .{ .Momentum, .AngularMomemtum } },
    // // .LinearDensity
    // .{
    //     .{ .Length, .Mass },
    // },
    // // .LuminousFlux
    // .{.{ .Time, .Energy }},
    // // .LuminousIntensity
    // .{
    //     // todo: add solid angle
    //     .{ .SolidAngle, .LuminousFlux },
    // },
    // // .MagneticFlux
    // .{
    //     // todo: Inverse here
    //     // .{InverseInductance, .Current}
    // },
    // // .Mass
    // .{
    //     .{ .Acceleration, .Force },
    //     .{ .Velocity, .Momemtum },
    //     .{ .LatentHeat, .Energy },
    // },
    // // .MeanLifetime = ,
    // // .Momemtum
    // .{
    //     // Add angular momentum
    //     .{ .Length, .AngularMomemtum },
    //     .{ .Velocity, .Energy },
    // },
    // // .MolarConcentration
    // .{
    //     .{ .Volume, .AmountOfSubstance },
    //     .{ .ChemicalPotential, .EnergyDensity },
    //     .{ .MolarEnergy, .EnergyDensity },
    // },
    // // .MolarEnergy
    // .{
    //     .{ .AmountOfSubstance, .Energy },
    //     .{ .MolarConcentration, .EnergyDensity },
    // },
    // // .MolarEntropy
    // .{
    //     .{ .AmountOfSubstance, .Entropy },
    // },
    // // .MolarHeatCapacity
    // .{
    //     .{ .AmountOfSubstance, .HeatCapacity },
    // },
    // // .MomentOfInertia
    // .{},
    // // .OpticalPower
    // .{},
    // // .Permeability
    // .{},
    // // .Permittivity
    // .{},
    // // .PotentialEnergy
    // .{},
    // // .Power
    // .{
    //     .{ .Frequency, .Energy },
    // },
    // // todo: should i descriminate between energy, work and heat?
    // // .Pressure
    // .{
    //     .{ .Volume, .Energy },
    //     .{ .Area, .Force },
    //     // todo: check that dynamic viscosity exists
    //     .{ .Time, .DynamicViscosity },
    //     .{ .Length, .SurfaceTension },
    //     // todo: strain, compressibility
    //     .{ .Compressibility, .Strain },
    // },
    // // .RadioActivity
    // .{},
    // // .RadiationDose
    // .{.{ .Mass, .Energy }},
    // // .Radiance
    // .{
    //     .{ .Area, .RadiantIntensity },
    //     .{ .SolidAngle, .Irradiance },
    // },
    // // .RadiantIntensity
    // .{.{ .SolidAngle, .Power }},
    // // .ReactionRate
    // .{},
    // // .Reluctance
    // .{
    //     .{ .MagneticFlux, .ElectricCurrent },
    // },
    // // .SpecificEnergy
    // .{
    //     .{ .Mass, .Energy },
    // },
    // // .SpecificHeatCapacity
    // .{
    //     .{ .Mass, .HeatCapacity },
    // },
    // // .SpecificVolume
    // .{ .Mass, .Volume },
    // // .Speed
    // .{},
    // // .Spin
    // .{.{ .Action, .AngularMomentum }},

    // // .Stress
    // .{
    //     .{ .Volume, .Energy },
    //     .{ .Area, .Force },
    //     .{ .Strain, .EnergyDensity },
    //     .{ .Time, .Viscosity },
    // },
    // // todo: check strain exists and what is it
    // // .Strain
    // .{
    //     .{ .YoungsModulous, .Stress },
    // },
    // // .SurfaceTension
    // .{
    //     .{ .Length, .Force },
    //     .{ .Area, .Energy },
    // },
    // // .Temperature
    // .{},
    // // .ThermalConductance
    // .{
    // // todo: Add temperature diff
    // .{ .TemperatureDiff, .Power }},
    // // .ThermalConductivity
    // .{
    //     .{ .Length, .ThermalConductance },
    // },
    // // .ThermalResistance
    // .{ .{ .Power, .TemperatureDiff }, .{ .Length, .ThermalResistivity } },
    // // .ThermalResistivity
    // .{},
    // // .Time
    // .{
    //     // Add impulse
    //     .{ .Force, .Impulse },
    //     .{ .Power, .Energy },
    //     .{ .Current, .EletricCharge },
    //     .{ .Acceleration, .Velocity },
    //     .{ .Frequency, .NumberOfCycles },
    // },
    // // .Viscosity
    // .{
    //     // TODO: this is a division actually
    //     // .{.InverseTime, .Stress},
    // },
    // // .Velocity
    // .{
    //     // Check displacement
    //     .{ .Time, .Length },
    //     .{ .Mass, .Momentum },
    //     .{ .Force, .Power },
    // },
    // // .Volume
    // .{
    //     .{ .Pressure, .Word },
    //     .{ .Density, .Mass },
    // },
    // // .VolumetricFlowRate
    // .{},
    // // .Wavelength
    // .{
    //     .{ .Frequency, .Speed },
    // },
    // // .Wavenumber
    // .{},
    // // .YoungsModulus
    // .{
    //     .{ .Strain, .Stress },
    // },
    // // .SpringConstant
    // .{.{ .Length, .Force }},
    // // todo: add solid angle
    // // .SolidAngle
    // .{},