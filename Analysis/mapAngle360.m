function mappedAngle = mapAngle360(angle)
    % Add 360 to negative angles
    angle(angle < 0) = angle(angle < 0) + 360;
    
    % Map angles greater than or equal to 360 to the range 0-360
    mappedAngle = mod(angle, 360);
end
