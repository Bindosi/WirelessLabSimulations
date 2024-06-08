function ind=dtct_edge2(sgnl_rx)
wnd_sz=8*16;
ab=abs(sgnl_rx).^2;
vv=vec2mat(ab,wnd_sz);
dd=diff(sum(vv.'));
plot(dd)
[~,ind]=max(dd);
ind=ind*wnd_sz;
end
