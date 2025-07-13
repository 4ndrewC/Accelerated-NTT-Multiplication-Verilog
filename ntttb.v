`timescale 1ns/1ns
`include "ntt.v"

module computetb;

    reg clk;
    reg [`WORD-1:0] cycles;
    reg rst;
    wire done;

    compute uut(clk, rst, cycles, done);

    always begin 
        clk = ~clk; #10;
    end

    always begin cycles = cycles + 1; #20; end


    initial begin
        $dumpfile("compute.vcd");
        $dumpvars(0, computetb);

        $display("\nSIMULATION BEGIN\n\n");
        
        cycles <= 0;
        clk <= 1'b1;
        rst = 1;
        #5;
        rst = 0;
        #1500;
        
        $finish;
    end;


endmodule
