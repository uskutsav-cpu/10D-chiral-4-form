def elimination(matrix,p):
    a=[list(map(lambda x:int(x)%p,row)) for row in matrix]
    columns=[];det=1
    for c in range(len(a[0])):
        r=len(columns)
        pivot=next((i for i in range(r,len(a)) if a[i][c]),None)
        if pivot is None:continue
        if pivot!=r:a[r],a[pivot]=a[pivot],a[r];det=-det
        v=a[r][c];det=det*v%p;inverse=pow(v,-1,p)
        for i in range(r+1,len(a)):
            scale=a[i][c]*inverse%p
            if scale:a[i]=[(x-scale*y)%p for x,y in zip(a[i],a[r])]
        columns.append(c)
        if len(columns)==len(a):break
    return columns,det if len(columns)==len(a) else 0


def orbit_rows(coordinates,p):
    lookup=dict(zip([t for t in TUPLES if 0 in t],coordinates))
    def component(indices):
        if len(set(indices))<5:return 0
        t=tuple(sorted(indices));s=parity(indices)
        if 0 in t:return s*lookup[t]%p
        complement=tuple(i for i in range(10) if i not in t)
        return -s*parity(t+complement)*lookup[complement]%p
    rows=[]
    for i,j in itertools.combinations(range(10),2):
        reverse=1 if i==0 else -1
        full=[]
        for t in TUPLES:
            total=0
            for k,index in enumerate(t):
                v=list(t)
                if index==i:v[k]=j;total+=component(v)
                elif index==j:v[k]=i;total+=reverse*component(v)
            full.append(total%p)
        row=full[:126]
        assert full==compact_form(row,p)
        rows.append(row)
    return rows


def validate_input(data):
    assert data['schema']==1
    assert data['convention']=={'metric':[-1]+[1]*9,'coordinate_tuples':[list(t) for t in TUPLES if 0 in t],'hodge':'output-first epsilon; epsilon_0123456789=+1; F=sum A_I(e_I+star e_I)'}
    assert len(data['graphs'])==len({g['id'] for g in data['graphs']})==81
    for item in data['graphs']:
        n=item['degree'];edges=item['graph']['edges']
        assert n==item['graph']['n'] and n>=4
        assert edges==sorted(edges) and len({(i,j) for i,j,m in edges})==len(edges)
        valences=[0]*n
        for i,j,m in edges:
            assert all(type(x)is int for x in (i,j,m)) and 0<=i<j<n and 0<m<=5
            valences[i]+=m;valences[j]+=m
        assert valences==[5]*n
    assert data['points']
    for point in data['points']:
        p=point['prime'];assert type(p)is int and 2<p<65536
        assert all(p%d for d in range(2,int(p**.5)+1))
        assert len(point['coordinates'])==126 and all(type(a)is int and 0<=a<p for a in point['coordinates'])
    for t in TUPLES:
        u,s=star_image(t);v,r=star_image(u)
        assert t==v and s*r==1


def main():
    if not __debug__:raise RuntimeError('Assertions must be enabled; do not use python -O')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();data=json.loads(args.input.read_text());validate_input(data)
    report={'schema':1,'input_sha256':hashlib.sha256(args.input.read_bytes()).hexdigest(),'implementation_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'mode':'full fresh values, Jacobians, Bareiss minors and orbit action from minimal inputs','points':[]}
    for point in data['points']:
        start=time.monotonic();p=point['prime'];a=point['coordinates'];engine=PointEngine(a,p)
        values=[];rows=[]
        for i,item in enumerate(data['graphs']):
            value,row=reverse_gradient(engine,item)
            assert sum(x*y for x,y in zip(a,row))%p==item['degree']*value%p
            values.append(value);rows.append(row)
            print(f'point prime={p} row={i+1}/81 id={item["id"]}',flush=True)
        columns,det=elimination(rows,p);assert len(columns)==81 and det!=0
        integer,_=bareiss([[row[c] for c in columns] for row in rows]);assert integer%p==det
        orbit=orbit_rows(a,p);oc,od=elimination(orbit,p);assert len(oc)==45 and od!=0
        orbit_integer,_=bareiss([[row[c] for c in oc] for row in orbit]);assert orbit_integer%p==od
        assert all(sum(x*y for x,y in zip(row,tangent))%p==0 for row in rows for tangent in orbit)
        report['points'].append({**point,'values':values,'jacobian':rows,'pivot_columns':columns,'integer_minor_determinant':str(integer),'determinant_mod_p':det,'rank':81,'orbit_rows':orbit,'orbit_pivot_columns':oc,'orbit_integer_minor_determinant':str(orbit_integer),'orbit_determinant_mod_p':od,'orbit_rank':45,'seconds':time.monotonic()-start,'arithmetic_checks':dict(engine.stats)})
    report['passed']=True
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'passed':True,'points':[{k:r[k] for k in ('prime','rank','determinant_mod_p','orbit_rank','orbit_determinant_mod_p','seconds')} for r in report['points']]},indent=2))
if __name__=='__main__':main()
