function manifest = buildRunManifest(projectDir,p,opts,exitflag,output)
%BUILDRUNMANIFEST Capture enough context to reproduce and audit a run.

manifest.schemaVersion = 1;
manifest.createdUTC = char(datetime('now','TimeZone','UTC', ...
    'Format','yyyy-MM-dd''T''HH:mm:ss.SSSXXX'));
manifest.matlabVersion = version;
manifest.toolboxes = ver;
manifest.parameters = p;
manifest.initialRandomSeed = p.randomSeed;
manifest.robustRandomSeed = p.robustSeed;
manifest.optimizerOptions = opts;
manifest.exitflag = exitflag;
manifest.solverOutput = output;
manifest.rngState = rng;

files = dir(fullfile(projectDir,'*.m'));
source = repmat(struct('name','','bytes',0,'modified','','sha256',''),numel(files),1);
for i = 1:numel(files)
    source(i).name = files(i).name;
    source(i).bytes = files(i).bytes;
    source(i).modified = files(i).date;
    source(i).sha256 = sha256File(fullfile(files(i).folder,files(i).name));
end
manifest.sourceFiles = source;
end

function hash = sha256File(path)
fid = fopen(path,'rb');
if fid<0
    hash = '';
    return
end
cleanup = onCleanup(@() fclose(fid));
bytes = fread(fid,Inf,'*uint8');
try
    md = java.security.MessageDigest.getInstance('SHA-256');
    md.update(typecast(bytes,'int8'));
    digest = typecast(md.digest(),'uint8');
    hash = lower(reshape(dec2hex(digest,2).',1,[]));
catch
    hash = 'unavailable';
end
end
