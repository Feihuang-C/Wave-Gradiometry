
function Anisomat2gmt(outfile,outpath,isbootstrap)
load(outfile)
periods=Aniso.period;
dM=1;
for pdi=1:length(periods)
    period=periods(pdi);
    st=Aniso.st;
    viso=Aniso.viso(:,pdi);
    stdv=Aniso.stdv(:,pdi);
    stdv2=Aniso.stdv2(:,pdi);
    Fai=Aniso.Fai(:,pdi);
    M=Aniso.M(:,pdi);
    if isbootstrap==1
        stdFai=Aniso.stdFai(:,pdi);
        stdM=Aniso.stdM(:,pdi);
    end

    Vnan = isnan(viso);
    st(Vnan,:)=[];
    viso(Vnan)=[];
    stdv(Vnan)=[];
    stdv2(Vnan)=[];


    Fai(Vnan)=[];
    M(Vnan)=[];
    
    if isbootstrap==1
        stdFai(Vnan)=[];
        stdM(Vnan)=[];
    end

    
    Vnan = find(isnan(Fai));
    st2=st;
    st2(Vnan,:)=[];
    stdv(Vnan)=[];    
    stdv2(Vnan)=[];    

    Fai(Vnan)=0;
    M(Vnan)=0;
    if isbootstrap==0
        fid=fopen([outpath,'V_prd_' num2str(period),'.txt'],'w');
        for i=1:length(M)
            fprintf(fid,'%f %f %f %f %f\n',st(i,1),st(i,2),viso(i),M(i),Fai(i));
        end
        fclose(fid);
        
          fid1=fopen([outpath,'err_prd_' num2str(period),'time.txt'],'w');
        for i=1:length(stdM)
            fprintf(fid1,'%f %f %f\n',st2(i,1),st2(i,2),stdv(i));
        end
        fclose(fid1);
        
        fid2=fopen([outpath,'err_prd_' num2str(period),'azm.txt'],'w');
        for i=1:length(stdM)
            fprintf(fid1,'%f %f %f\n',st2(i,1),st2(i,2),stdv2(i));
        end
        fclose(fid2);
        
    else
        fid=fopen([outpath,'V_prd_' num2str(period),'.txt'],'w');
        for i=1:length(M)
            if stdM(i)<dM
                fprintf(fid,'%f %f %f %f %f\n',st(i,1),st(i,2),viso(i),M(i),Fai(i));
            else
                fprintf(fid,'%f %f %f %f %f\n',st(i,1),st(i,2),viso(i),0,0);
            end
        end
        fclose(fid);
        stdFai(Vnan)=[];
        stdM(Vnan)=[];
        fid1=fopen([outpath,'err_prd_' num2str(period),'time.txt'],'w');
        for i=1:length(stdM)
            fprintf(fid1,'%f %f %f %f %f\n',st2(i,1),st2(i,2),stdv(i),stdM(i),stdFai(i));
        end
        fclose(fid1);
        
        fid2=fopen([outpath,'err_prd_' num2str(period),'azm.txt'],'w');
        for i=1:length(stdM)
            fprintf(fid1,'%f %f %f %f %f\n',st2(i,1),st2(i,2),stdv2(i),stdM(i),stdFai(i));
        end
        fclose(fid2);
    
        
        
    end
end
