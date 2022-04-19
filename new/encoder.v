`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 14.04.2022 14:37:30
// Design Name: 
// Module Name: encoder
// Project Name: 

module encoder(
    input wire clk,
    input wire reset,
    input wire doneThreshold,
    input wire signed [11 : 0] threshold_sig,
    input wire [2047 : 0] prbs_in, 
    input wire [11 : 0] id,
    
//    output reg we, // RAM signals
//    output reg [8 : 0] addr,
//    output reg signed [16 : 0] y_out,
    output reg done
    );
    
parameter M = 512, LENGTH = 256;

wire [2 : 0] index;
assign index = id[10 : 8];

reg signed [16 : 0] y [0 : 511]; //Final RAM output
reg [7 : 0] cnt;
integer j;

always @ (posedge clk) begin
    if (reset == 1) begin
        for(j = 0; j < M; j = j + 1) begin
            y[j] <= 17'b0;
        end
        done <= 1'b0;  
        cnt <= 8'b0;
    end
    
    else begin
        if(doneThreshold == 1) begin
            for(j = 0; j < 64; j = j + 1) begin
                y[j + index] <= prbs_in[{index, 8'b0} + (j >> 3) + cnt] ? (y[j + index] + threshold_sig) : (y[j + index] - threshold_sig);
            end
            cnt <= cnt + 1;
        end
//        done <= 1'b1;
    end
end
endmodule
