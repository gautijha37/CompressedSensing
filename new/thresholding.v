`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 16.03.2022 10:18:33
// Design Name: 
// Module Name: thresholding
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

module thresholding(
    input wire clk, reset, enable,
    input wire signed [11 : 0] sig,
    
    output reg signed [11 : 0] y,
    output reg done,
    output reg [11 : 0] id
    );
    
parameter N = 2048, W = 13, threshold = 84, MEAN = 958;
localparam IDLE = 1'b0,
           WAIT = 2'b1,
           WORKING = 1'b1;
           
reg signed [11 : 0] d [0 : W - 1]; // Shift regs to hold W length of values
reg signed [11 : 0] th;
reg signed [16 : 0] sum;
reg signed [11 : 0] curr_avg;

wire signed [11 : 0] din_curr_avg;
wire signed [16 : 0] din_sum;

reg state;
wire signed [11 : 0] mean;
wire signed [11 : 0] sig_diff;

assign din_curr_avg = din_sum/W;
assign mean = MEAN;
assign sig_diff = sig - mean;
assign din_sum = (sum + sig_diff - d[W - 1]);

integer i, j;
//assign done = (state == WORKING);
//always @ posedge clk
//begin
//    if(reset) begin
//        d <= 0;
//    end
    
//    else begin
//        d[1 : W - 1] <= d[0 : W - 2];
//        d[0] <= sig;
//    end
//end


always @ (posedge clk)
begin
    if(reset) begin
        for(i = 0; i < W; i = i + 1) begin
            d[i] <= 12'b0;
        end
        sum <= 17'b0;
        curr_avg <= 12'b0;
        th <= threshold;
        state <= IDLE;
        done <= 1'b0;
        id <= 12'b0;
    end
    
    else begin
        case(state)
        IDLE: begin
            if(enable == 1'b1) begin
                state <= WORKING;
                id <= 12'b111111111111;
            end
            done <= 1'b0;
            
        end
        
//        WAIT : begin
//            state <= WORKING;
//        end
        
        WORKING: begin
            sum <= din_sum;
            curr_avg <= din_curr_avg;
//            for(j = 0;
//            d[1 : W - 1] <= d[0 : W - 2];
            for(j = 1; j < W; j = j + 1) begin
                d[j] <= d[j - 1];
            end
            d[0] <= sig_diff;
            
            if(sig_diff >= th + curr_avg) begin
                y <= sig_diff;
            end
            else begin
                y <= 12'b0;
            end
            
            done <= 1'b1;
            id = id + 1;
            
            if(id == 2047)
                state <= IDLE;
//            if(id == 
//            if(id == 0)
        end
        endcase
    end
end
endmodule
