`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 09.02.2022 16:27:18
// Design Name: 
// Module Name: tb_prbs
// Project Name: 
// Target Devices: 
// Tool Versions: 
// Description: 
// 
// Dependencies: 
// 
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// 
//////////////////////////////////////////////////////////////////////////////////


module tb_prbs;
//Inputs
reg clk, reset;

//Outputs
wire [2047 : 0] prbs;

//Module instantiation
prbs t1(.clk(clk), .reset(reset), 
        .prbs(prbs));
        
//Stimulus
initial
begin
    clk = 1'b0;
    repeat(5000)
    #5 clk = ~clk;
    $stop;
end

initial
begin
    #7 reset <= 1'b1;
    #10 reset <= 1'b0; 
end       

parameter outfile = "prbs.txt";
integer fid, j;

initial begin
    fid = $fopen (outfile, "w");
    j = 0;
end

always @ (posedge clk) begin
    if(reset == 1'b0) begin
        if(j < 2048) begin
            $fdisplay(fid, prbs[j]);
            j = j + 1;
        end
    end
end
/*
//Read outputs from text file
parameter ECG_FILE = "ecg.mif";
reg [31:0] outa1 [0 : 11];

integer fid_out1, j1;
initial begin
    fid_out1 = $fopen (ECG_FILE, "r");
    j1 = 0;
    while(! $feof(fid_out1)) begin
        $fscanf(fid_out1, "%d\n", outa1[j1]);
        j1 = j1 + 1;
    end    
    $fclose(fid_out1);
end
*/


endmodule
