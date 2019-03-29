function writeOff(fname,X,T)

fp = fopen(fname,'w');
fwrite(fp,sprintf('OFF\n'));
fwrite(fp,sprintf('%d %d 0\n',size(X,1),size(T,1)));
dlmwrite(fname,X,'-append','delimiter',' ');
dlmwrite(fname,[repmat(3,size(T,1),1) T-1],'-append','delimiter',' ');
fclose(fp);

end