function Talb=IterTalb(T,So)
% assumes T has been generated using ebm and that ebmAlbedo has the
% same parameters (D,albedo, A, B) as ebm.

diffoff=1;
Talb=T;
counter=0;
while (diffoff > 1e-5)
  [Talb_new,xm]=ebmAlbedo(Talb,So);
  diffoff=sum(abs(Talb_new-Talb));
  counter=counter+1;
  Talb=Talb_new;
  if (counter >50) break; end
end
disp(['solar constant ',num2str(So),' counter ',num2str(counter),' difference ',num2str(diffoff*1e5),'E-5'])

return
  
