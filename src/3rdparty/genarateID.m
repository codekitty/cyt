%Unique id generator 
function uniqueID = genarateID
    rng('shuffle')
    symbolsLetters = ['a':'z' 'A':'Z'];
    symbolsNumbers = ['0':'9'];
    nums = randi(numel(symbolsNumbers),[1 4]);
    lets = randi(numel(symbolsLetters),[1 4]);
    str = symbolsLetters (lets);
    nums=num2str(nums);
    nums = nums(nums ~= ' ');
    uniqueID= [str,nums];
end