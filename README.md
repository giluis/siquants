Siquants: SI unit conversion and math library

Siquants (seequants) provides a generic struct Q which represents a physical quantity like length, eletrical current or mass.

## Features:
- Compile time validated addition of the same (or identical) quantities, 
  - ✅ Length + Length -> Length
  - ✅ Work   + Heat -> Energy   
  - ⛔️ Length + Mass --> compile error
- Compile time validated multiplication of sensible SI quantities
  - ✅ Length * Length = Area       
  - ✅ Force  / Area   = Pressure   
  - ⛔️ Electrical Capactity / Refractive Index --> compile error  
- Quantities are represented as just an `f64`, no extra space needed
- You can decide what is the base SI order of magnitude you want the quantity to be represented
  - This allows quantities to be centered around the ordered of magnitude of your preference, allowing you to squeeze as much precision as possible for the values you are trying to represent
  - User should track precision loss themselves

## Instalation
TODO

## Quick start:
```zig  

const print = std.debug.print;
const Q = @import("siquants").Quantity;

// Quantity is a type defined like:
//
// fn Quantity (comptime quantkind: QuantityKind, comptime base: SIOrder) type {...}
//
// where
// - QuantityKind is one of the many Physical Quantity kinds in the SI system,
// - base is something like .Milli or .Tera, which will maximize the
//    precision of the floating point value around that order of of magnitude

fn main() void {
  // initialize a Length as 4.0 Millis.
  const len1 = Q(.Length, .Milli).init(4.0);

  print("len1 as kilometers {}", .{len1.as(.Kilo).value});

  // 24 km stored as milimiters. Might lose some precision
  const len2 = Q(.Length, .Milli).from(.Kilo, 24.0);

  const area = len1.mult(len2);

  // Convert to meters for printing
  print("A rectangle of 4mm * 24km has an area of {} m^2", .{area.as(.Unit)});
}

```

