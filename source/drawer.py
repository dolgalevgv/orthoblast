import os
import copy
import pygame as pg
from pygame.locals import *

from . import misc

def drawResult(spec, prot, ev, cov, result, outname):
    def drawRect(x1, x2, y, color):
        pg.draw.rect(screen, misc.COLOR[color], Rect(x1, y, x2-x1, 14))
    def prepHit(hit):
        chit = copy.deepcopy(hit)
        gmin = 100000000000
        gmax = 0
        for h in chit.hsps:
            h.qstart /= chit.query.qlen / 500
            h.qend /= chit.query.qlen / 500
            h.qstart = round(h.qstart)
            h.qend = round(h.qend)
            if h.gstart < gmin:
                gmin = h.gstart
            if h.gend > gmax:
                gmax = h.gend
        for h in chit.hsps:
            h.gstart -= gmin
            h.gstart /= (gmax - gmin) / 500
            h.gend -= gmin
            h.gend /= (gmax - gmin) / 500
            h.gstart = round(h.gstart)
            h.gend = round(h.gend)
        if chit.domains != []:
            for d in chit.domains:
                d.start /= chit.query.qlen / 500
                d.end /= chit.query.qlen / 500
        for i in range(1, len(chit.hsps)):
            if hit.hsps[i].qstart - hit.hsps[i-1].qend == 1:
                chit.hsps[i].qstart = chit.hsps[i-1].qend - 2
        return chit
    def prepQuery(query):
        q = copy.deepcopy(query)
        for d in q.domains:
            d.start /= q.qlen / 500
            d.end /= q.qlen / 500
        return q
    def drawTitle():
        title = font.render(f'HUMAN {prot.upper()} VS {spec.upper()} GENOME (E={ev if ev else 10}, COV={cov if cov else "CUSTOM"})', True, misc.COLOR['black'])
        screen.blit(title, (10, 7))
    def drawHeatMap():
        for i in range(1, 9):
            pg.draw.rect(screen, misc.COLOR[f'hm{i}'], Rect(370+i*24, 6, 24, 10))
            score = font.render(str(30*i), True, misc.COLOR['black'])
            screen.blit(score, (382+i*24-score.get_width()/2, 18))
    def drawDomain(x1, x2, y, name):
        drawRect(x1, x2, y, 'e')
        dname = font.render(name, True, misc.COLOR['black'])
        coord = (x1 + (x2-x1)/2) - dname.get_width()/2
        screen.blit(dname, (coord, y+2))
    def drawQuery(query):
        query = prepQuery(query)
        drawRect(50, 550, 32, 'c')
        qname = font.render('QUERY', True, misc.COLOR['black'])
        qsize = font.render(str(query.qlen)+'bp', True, misc.COLOR['black'])
        screen.blit(qname, (25-qname.get_width()/2, 34))
        screen.blit(qsize, (575-qsize.get_width()/2, 34))
        for d in query.domains:
            drawDomain(50+d.start, 50+d.end, 32, d.name)
    def drawHit(y, hit, flag):
        hit = prepHit(hit)
        for hsp in hit.hsps:
            if hsp.score > 210:
                color = 'hm8'
            elif hsp.score > 180:
                color = 'hm7'
            elif hsp.score > 150:
                color = 'hm6'
            elif hsp.score > 120:
                color = 'hm5'
            elif hsp.score > 90:
                color = 'hm4'
            elif hsp.score > 60:
                color = 'hm3'
            elif hsp.score > 30:
                color = 'hm2'
            else:
                color = 'hm1'
            drawRect(50+hsp.qstart, 50+hsp.qend, y, color)
            drawRect(50+hsp.gstart, 50+hsp.gend, y+20, color)
        nhit = font.render('HIT', True, misc.COLOR['black'])
        npos = font.render('POS', True, misc.COLOR['black'])
        nseq = font.render('SEQ', True, misc.COLOR['black'])
        screen.blit(nhit, (25-nhit.get_width()/2, y+2))
        screen.blit(nseq, (575-nseq.get_width()/2, y+2))
        screen.blit(npos, (575-npos.get_width()/2, y+22))
        if flag:
            y += 40
            for d in hit.domains:
                drawDomain(50+d.start, 50+d.end, y, d.name)
            ndom = font.render('DOM', True, misc.COLOR['black'])
            screen.blit(ndom, (575-ndom.get_width()/2, y+2))
    pg.init()
    length = 60*len([hit for hit in result.hits if hit.seq != '']) + 40*len([hit for hit in result.hits if hit.seq == ''])
    screen = pg.display.set_mode((600, 60+length))
    screen.fill((245,245,245))
    font = pg.font.SysFont('Arial Bold', 15)
    drawTitle()
    drawHeatMap()
    drawQuery(result.query)
    y = 60
    for i in range(len(result.hits)):
        if result.hits[i].seq != '':
            drawHit(y, result.hits[i], True)
            y += 60
        else:
            drawHit(y, result.hits[i], False)
            y += 40

    if not os.path.isdir(os.path.join(os.getcwd(), 'results', f'{outname if outname else ""}', 'graphs')):
        os.mkdir(os.path.join(os.getcwd(), 'results', f'{outname if outname else ""}', 'graphs'))
    pg.image.save(screen, os.path.join(os.getcwd(), 'results', f'{outname if outname else ""}', 'graphs', f'{spec}_{prot}.png'))
    pg.quit()
