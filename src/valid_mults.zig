const QuantityKind = @import("quantkind.zig").QuantityKind;
const std = @import("std");
const MultPair: type = [2]QuantityKind;

pub fn a() void {
    const mult_fields = @typeInfo(VALID_MULTS).@"struct".fields;
    const qk_fields = @typeInfo(QuantityKind).@"enum".fields;
    // comptime std.debug.assert(mult_fields.len == qk_fields.len);
    @setEvalBranchQuota(20000);

    std.debug.print("\n\n ===== Quantity Kind is lacking... ", .{});
    inline for (mult_fields) |m| {
        inline for(qk_fields) |q| {
            var lower_m = [_]u8{undefined} ** 100;
            var lower_q = [_]u8{undefined} ** 100;
            _ = std.ascii.lowerString(&lower_m, m.name);
            _ = std.ascii.lowerString(&lower_q, q.name);
            if (std.mem.eql(u8, &lower_m, &lower_q)) {
                break;
            }
        } else {
            std.debug.print("\n- {s}", .{m.name});

        }
    }
    std.debug.print("\n\n ===== Valid Mults is lacking... ", .{});

    inline for(qk_fields) |q| {
        inline for (mult_fields) |m| {
            var lower_m = [_]u8{undefined} ** 100;
            var lower_q = [_]u8{undefined} ** 100;
            _ = std.ascii.lowerString(&lower_m, m.name);
            _ = std.ascii.lowerString(&lower_q, q.name);
            if (std.mem.eql(u8, &lower_m, &lower_q)) {
                break;
            }
        } else {
            std.debug.print("\n- {s}", .{q.name});
        }
    }
}

const Operation = enum {
    Div, 
    Mul, 
}; 

const OperationGraph = struct {
    const num_quantity_kinds = @typeInfo(QuantityKind).@"enum".fields.len;
    quantity_kinds:  [num_quantity_kinds] QuantityKind, 
    edges: [num_quantity_kinds] [10] Operation,
};

pub const VALID_MULTS = ValidMultsStruct{}; 

// TODO: check this is correct
const ValidMultsStruct = struct {
    AbsorbedDoseRate: [2]MultPair = .{ .{ .Mass, .Power }, .{ .HalfLife, .RadiationDose } },
    Action: [1]MultPair = .{.{ .Frequency, .Energy }},
    AmountOfSubstance: [4]MultPair = .{
        .{ .ChemicalPotential, .Energy },
        .{ .MolarEnergy, .Energy },
        .{ .MolarEntropy, .Entropy },
        .{ .MolarHeatCapacity, .HeatCapacity },
    },
    AngularAcceleration: [0]MultPair = .{},
    // todo: add mults for angular momentum
    AngularMomemtum: [0]MultPair = .{},
    Acceleration: [3]MultPair = .{
        .{ .Mass, .Force },
        .{ .Time, .Velocity },
        // Add Energy per mass
        .{ .Length, .EnergyPerMass },
    },
    Area: [5]MultPair = .{
        .{ .Length, .Volume },
        .{ .Pressure, .Force },
        .{ .Illuminance, .LuminousFlux },
        .{ .MagneticField, .MagneticFlux },
        // todo(this could be work)
        .{ .SurfaceTension, .Work },
    },
    AreaDensity: [1]MultPair = .{.{ .Area, .Mass }},
    Capacitance: [1]MultPair = .{
        .{ .ElectricPotential, .ElectricCharge },
        // todo: consider Voltage^2 here. Also, does voltage^2 exist here?
        // .{.ElectricPotentialSquared, .ElectricCharge},
    },
    Compressibility: [0] MultPair = .{},
    // todo: add Vector
    MagneticField: [1]MultPair = .{
        .{ .Area, .MagneticFlux },
    },
    CatalyticActivityConcentration: [0]MultPair = .{},
    ChemicalPotential: [1]MultPair = .{
        .{ .MolarConcentration, .EnergyDensity },
    },
    VolumeDensity: [1]MultPair = .{
        .{ .Volume, .Mass },
    },
    // todo: add density
    Density: [0] MultPair = .{},
    DynamicViscosity: [0] MultPair = .{},
    DoseEquivalent: [1]MultPair = .{.{ .Mass, .Energy }},
    ElectricCharge: [3]MultPair = .{
        // todo: consider voltage herecross product of velocity and idk what else, need help here
        .{ .ElectricPotential, .Energy },
        // todo: add electric field
        .{ .ElectricField, .Force },
        // todo: this is the cross product of velocity and idk what else, need help herecross product of velocity and idk what else, need help here
        .{ .MagneticField, .Force },
    },
    ElectricChargeDensity: [0]MultPair = .{},
    ElectricalConductance: [1]MultPair = .{
        .{ .ElectricPotential, .ElectricCurrent },
    },
    ElectricalConductivity: [0]MultPair = .{},
    ElectricCurrent: [3]MultPair = .{
        .{ .Time, .ElectricCharge },
        // todo: consider voltage here
        .{ .ElectricalResistance, .ElectricPotential },
        // todo: consider voltage here
        .{ .ElectricPotential, .Power },
    },
    ElectricField: [2]MultPair = .{
        .{ .ElectricCharge, .Force },
        // todo: consider electric potential
        .{ .Length, .ElectricPotential },
    },
    // todo: consider voltage here
    ElectricPotential: [4]MultPair = .{
        .{ .ElectricCurrent, .Power },
        .{ .ElectricCharge, .Energy },
        .{ .ElectricalConductance, .ElectricCurrent },
        .{ .Time, .MagneticFlux },
    },
    ElectricalResistance: [2]MultPair = .{
        .{ .ElectricCurrent, .ElectricPotential },
        .{ .Length, .ElectricalResistivity },
        // todo: add a name for this
        // .{ .ElectricCurrentSquared, .ElectricPotential },
    },
    ElectricalResistivity: [0]MultPair = .{},
    Energy: [2]MultPair = .{ .{ .Frequency, .Power }, .{ .Time, .Action } },
    // todo: energy per mass mults
    EnergyPerMass: [0]MultPair = .{},
    EnergyDensity: [1]MultPair = .{.{ .Volume, .Energy }},
    Entropy: [1]MultPair = .{
    // todo: should I have heat
    .{ .Temperature, .Heat }},
    Frequency: [1]MultPair = .{
        // TODO: add physical constants
        // .{ .PlancksConstant, .Energy },
        // todo: add number of cycles
        .{ .Time, .NumberOfCycles },
    },
    // todo: Vector
    Force: [3]MultPair = .{
        .{ .Time, .Impulse },
        // todo: Should I separate length and distance
        // Applicable for Torque as well
        .{ .Length, .Energy },
        .{ .Velocity, .Power },
    },
    HalfLife: [1]MultPair = .{.{ .AbsorbedDoseRate, .RadiationDose }},
    Heat: [0]MultPair = .{},
    // heat = ,=
    HeatCapacity: [1]MultPair = .{
    // todo: add heat?
    .{ .Temperature, .Heat }},
    HeatFluxDensity: [0]MultPair = .{},
    Illuminance: [1]MultPair = .{.{ .Area, .LuminousFlux }},
    Impedance: [1]MultPair = .{.{ .ElectricCurrent, .ElectricPotential }},
    // todo: impulse mults
    Impulse: [0]MultPair = .{},
    Inductance: [1]MultPair = .{
        .{ .ElectricCurrent, .MagneticFlux },
        // todo: not sure if should add current squared
        // .{ .ElectricCurrentSquared, .MagneticFlux },
    },
    Irradiance: [1]MultPair = .{.{ .Area, .Power }},
    Intensity: [1]MultPair = .{.{ .Area, .Power }},
    KineticEnergy: [0]MultPair = .{},
    // kineticEnergy = ,=
    // todo: should I have latent heat?
    LatentHeat: [0]MultPair = .{},
    Length: [4]MultPair = .{ .{ .Length, .Area }, .{ .Area, .Volume }, .{ .Force, .Energy }, .{ .Momemtum, .AngularMomemtum } },
    LinearDensity: [1]MultPair = .{
        .{ .Length, .Mass },
    },
    LuminousFlux: [1]MultPair = .{.{ .Time, .Energy }},
    LuminousIntensity: [1]MultPair = .{
        // todo: add solid angle
        .{ .SolidAngle, .LuminousFlux },
    },
    MagneticFlux: [0]MultPair = .{
        // todo: Inverse here
        // .{.InverseInductance, .Current}
    },
    Mass: [3]MultPair = .{
        .{ .Acceleration, .Force },
        .{ .Velocity, .Momemtum },
        .{ .LatentHeat, .Energy },
    },
    // meanLifetime = ,=
    Momemtum: [2]MultPair = .{
        // Add angular momentum
        .{ .Length, .AngularMomemtum },
        .{ .Velocity, .Energy },
    },
    MolarConcentration: [3]MultPair = .{
        .{ .Volume, .AmountOfSubstance },
        .{ .ChemicalPotential, .EnergyDensity },
        .{ .MolarEnergy, .EnergyDensity },
    },
    MolarEnergy: [2]MultPair = .{
        .{ .AmountOfSubstance, .Energy },
        .{ .MolarConcentration, .EnergyDensity },
    },
    MolarEntropy: [1]MultPair = .{
        .{ .AmountOfSubstance, .Entropy },
    },
    MolarHeatCapacity: [1]MultPair = .{
        .{ .AmountOfSubstance, .HeatCapacity },
    },
    // todo: momentofinertial mults
    MomentOfInertia: [0]MultPair = .{},
    // todo: numberofcycles mults
    NumberOfCycles: [0]MultPair = .{},
    OpticalPower: [0]MultPair = .{},
    Permeability: [0]MultPair = .{},
    Permittivity: [0]MultPair = .{},
    // TODO(luis.gil): add this
    // potentialEnergy:  [_]MultPair=.{},
    Power: [1]MultPair = .{
        .{ .Frequency, .Energy },
    },
    // todo: should i descriminate between energy, work and heat?
    Pressure: [5]MultPair = .{
        .{ .Volume, .Energy },
        .{ .Area, .Force },
        // todo: check that dynamic viscosity exists
        .{ .Time, .DynamicViscosity },
        .{ .Length, .SurfaceTension },
        // todo: strain, compressibility
        .{ .Compressibility, .Strain },
    },

    RadioActivity: [0]MultPair = .{},
    RadiationDose: [1]MultPair = .{.{ .Mass, .Energy }},
    Radiance: [2]MultPair = .{
        .{ .Area, .RadiantIntensity },
        .{ .SolidAngle, .Irradiance },
    },
    RadiantIntensity: [1]MultPair = .{.{ .SolidAngle, .Power }},
    ReactionRate: [0]MultPair = .{},
    Reluctance: [1]MultPair = .{
        .{ .MagneticFlux, .ElectricCurrent },
    },
    SpecificEnergy: [1]MultPair = .{
        .{ .Mass, .Energy },
    },
    SpecificHeatCapacity: [1]MultPair = .{
        .{ .Mass, .HeatCapacity },
    },
    SpecificVolume: [1]MultPair = .{.{ .Mass, .Volume }},
    Speed: [0]MultPair = .{},
    Spin: [1]MultPair = .{.{ .Action, .AngularMomemtum }},
    Stress: [4]MultPair = .{
        .{ .Volume, .Energy },
        .{ .Area, .Force },
        .{ .Strain, .EnergyDensity },
        .{ .Time, .Viscosity },
    },
    // todo: check strain exists and what is it
    Strain: [1]MultPair = .{
        .{ .YoungsModulous, .Stress },
    },
    SurfaceTension: [2]MultPair = .{
        .{ .Length, .Force },
        .{ .Area, .Energy },
    },
    Temperature: [0]MultPair = .{},
    // todo: temperatureDiff mults
    TemperatureDiff: [0]MultPair = .{},
    ThermalConductance: [1]MultPair = .{
    // todo: Add temperature diff
    .{ .TemperatureDiff, .Power }},
    ThermalConductivity: [1]MultPair = .{
        .{ .Length, .ThermalConductance },
    },
    ThermalResistance: [2]MultPair = .{ .{ .Power, .TemperatureDiff }, .{ .Length, .ThermalResistivity } },
    ThermalResistivity: [0]MultPair = .{},
    Time: [5]MultPair = .{
        // Add impulse
        .{ .Force, .Impulse },
        .{ .Power, .Energy },
        .{ .ElectricCurrent, .ElectricCharge },
        .{ .Acceleration, .Velocity },
        .{ .Frequency, .NumberOfCycles },
    },
    Viscosity: [0]MultPair = .{
        // TODO: this is a division actually
        // .{.InverseTime, .Stress},
    },
    Velocity: [3]MultPair = .{
        // Check displacement
        .{ .Time, .Length },
        .{ .Mass, .Momemtum },
        .{ .Force, .Power },
    },
    Volume: [2]MultPair = .{
        .{ .Pressure, .Energy },
        .{ .Density, .Mass },
    },
    VolumetricFlowRate: [0]MultPair = .{},
    Wavelength: [1]MultPair = .{
        .{ .Frequency, .Speed },
    },
    Wavenumber: [0]MultPair = .{},
    Work: [0]MultPair = .{},
    YoungsModulous: [1]MultPair = .{
        .{ .Strain, .Stress },
    },
    SpringConstant: [1]MultPair = .{.{ .Length, .Force }},
    // todo: add solid angle
    SolidAngle: [0]MultPair = .{},
    // todo: make this unreachable
    UNKNOWN: [0]MultPair = .{},
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
// // .Acceleration
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
// .{.{ .Action, .AngularMomemtum }},

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
//     .{ .Current, .ElectricCharge },
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
