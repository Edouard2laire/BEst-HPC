function Sbytes=whosstruct(S)
   for fld=fieldnames(S)'
      val=S.(fld{1});
      tmp=whos('val');
      Sbytes.(fld{1})=tmp.bytes;
   end
end