clear,
% % Message Source and Bits
% Message = ['We hold these truths 2 be self-evident,' ...
%     ' that all men are created equal, that they are ' ...
%     'endowed by their Creator with certain unalienable ' ...
%     'Rights, that among these are Life, Liberty and ' ...
%     'the pursuit of Happiness'];
% 
% MessageBinChar = dec2bin(double(char(Message)));
% [Sz1,Sz2]      = size(MessageBinChar);
% 
% MessageBinStrn = reshape(MessageBinChar,1,Sz1*Sz2);
% MessageBits    = zeros(1,Sz1*Sz2);
% 
% for b=1:Sz1*Sz2
%     MessageBits(b) = str2double(MessageBinStrn(b));
% end

MessageBits = char2bin(['We hold these truths 2 be self-evident,' ...
     ' that all men are created equal, that they are ' ...
     'endowed by their Creator with certain unalienable ' ...
     'Rights, that among these are Life, Liberty and ' ...
     'the pursuit of Happiness']);

MessageBits = [MessageBits zeros(1,24)];
% Define your sequence of bits
bits = MessageBits; % Example bit sequence for "hello"
bits = bin2char(bits)

function origchar = bin2char(binvec)
   origchar = char(reshape(bin2dec(reshape(char(binvec + '0'),8,[]).'),1,[]));
end


function binvec = char2bin(charinput)
    binvec = reshape((dec2bin(charinput,8) - '0').',1,[]);
end