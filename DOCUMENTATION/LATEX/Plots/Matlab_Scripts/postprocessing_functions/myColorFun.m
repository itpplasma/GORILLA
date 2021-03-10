function [color] = myColorFun(color_in)

    switch color_in
        case 'green'
        	color = [18/255,173/255,42/255];   
        case 'grey'
            color = [0.7,0.7,0.7];
        case 'lightblue'
            color = [0,178,255]/255;
        case 'darkgrey'
            color = [0.3,0.3,0.3];
        case 'gold'
            color = [249,166,2]/255;
    end
    
end

