function divfunc(a, b)
    -- Function to calculate the mean of a table
    local function mean(t)
        local sum = 0
        local count = 0

        for k,v in pairs(t) do
            if type(v) == 'number' then
                sum = sum + v
                count = count + 1
            end
        end

        return (sum / count)
    end

    -- Check if a or b is a table and calculate the mean if necessary
    if type(a) == 'table' then
        a = mean(a)
    end
    if type(b) == 'table' then
        b = mean(b)
    end

    -- Original function
    if a == nil or b == nil then
        return ""
    elseif a == 0 then
        return "0.00"
    else
        return string.format("%d", a / b)
    end
end