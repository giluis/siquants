const std = @import("std");
const print = std.debug.print;

/// Will all non unique elements while using the arrayList itself as a HashSet
pub fn removeNonUnique(comptime T: type, list: *std.ArrayList(T), _: fn (*T) usize) void {
    const ptr = list.items.ptr;
    const cap = list.capacity;
    print("len: {}, {any}  ", .{ list.items.len, ptr[0 .. cap + 1000] });
    // print("len:  {any}  ", .{list.items[list.items.len]});
    print("\n {}", .{ptr[7]});
}

pub fn require_struct_field(Ty: type, field_name: []const u8, required_type: type) void {

    for (@typeInfo(Ty).@"struct".fields) |field| {
        if (std.mem.eql(u8, field.name, field_name) and field.type == f64) {
            break;
        }
        if (std.mem.eql(u8, field_name, @typeName(required_type))) {
            break;
        }
    } else {
        @compileError("please provide a \"" ++ field_name ++ "\" field of type f64");
    }
}

pub fn fieldIndex(comptime MyStruct: type, comptime field_name: []const u8) usize {
    return switch (@typeInfo(MyStruct)) {
        .@"struct" => |struct_info| inline for (struct_info.fields, 0..) |field, idx| {
            if (std.mem.eql(u8, field.name, field_name)) return idx;
        } else @compileError("field does not exist in struct"),
        else => @compileError("Please provide a struct"),
    };
}

pub fn OneWayPairIterator(comptime T: type) type {
    return struct {
        const Self = @This();
        data: []const T,
        i: usize,
        j: usize,

        pub fn init(data: []const T) Self {
            return Self{
                .i = 0,
                .j = 0,
                .data = data,
            };
        }

        pub fn next(self: *Self) ?[2]*const T {
            if (self.j == self.data.len - 1) {
                if (self.i == self.data.len - 2) {
                    return null;
                }
                self.i += 1;
                self.j = self.i + 1;
            } else {
                self.j += 1;
            }
            return .{ &self.data[self.i], &self.data[self.j] };
        }
    };
}

pub fn PairIterator(comptime T: type) type {
    return struct {
        const Self = @This();
        data: []const T,
        i: usize,
        j: usize,
        invert: bool,

        pub fn init(data: []const T) Self {
            return Self{
                .i = 0,
                .j = 0,
                .invert = false,
                .data = data,
            };
        }

        pub fn next(self: *Self) ?[2]*const T {
            if (self.invert) {
                self.invert = false;
                return .{ &self.data[self.j], &self.data[self.i] };
            }

            self.invert = true;

            if (self.j == self.data.len - 1) {
                if (self.i == self.data.len - 2) {
                    return null;
                }
                self.i += 1;
                self.j = self.i + 1;
            } else {
                self.j += 1;
            }
            return .{ &self.data[self.i], &self.data[self.j] };
        }
    };
}