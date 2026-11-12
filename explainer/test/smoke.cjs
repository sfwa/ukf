/*
Open the built site in a real browser and check that the article works:
every figure draws, nothing reports an error, every link in the text can be
pressed, and every number in the text is the number the estimator gives.

Usage:
    node explainer/test/smoke.cjs [SITE_DIRECTORY]

Needs the "playwright" package and its Chromium (npx playwright install
chromium). The site directory defaults to explainer/_site.
*/

const fs = require('node:fs');
const http = require('node:http');
const path = require('node:path');
const {chromium} = require('playwright');

const site = path.resolve(process.argv[2] ?? path.join(__dirname, '..', '_site'));
const TYPES = {
    '.html': 'text/html', '.js': 'text/javascript', '.css': 'text/css', '.json': 'application/json',
    '.wasm': 'application/wasm', '.pack': 'application/octet-stream',
};

/* A static file server for the site folder. */
function serve() {
    const server = http.createServer((request, response) => {
        const name = decodeURIComponent(new URL(request.url, 'http://localhost').pathname);
        const file = path.join(site, name.endsWith('/') ? `${name}index.html` : name);
        if (!file.startsWith(site) || !fs.existsSync(file) || fs.statSync(file).isDirectory()) {
            response.writeHead(404).end();
            return;
        }
        response.writeHead(200, {'Content-Type': TYPES[path.extname(file)] ?? 'application/octet-stream'});
        fs.createReadStream(file).pipe(response);
    });
    return new Promise(resolve => server.listen(0, '127.0.0.1', () => resolve(server)));
}

async function main() {
    const server = await serve();
    const url = `http://127.0.0.1:${server.address().port}/`;
    const browser = await chromium.launch();
    const problems = [];
    try {
        const page = await browser.newPage({viewport: {width: 1200, height: 900}, reducedMotion: 'reduce'});
        page.on('console', message => {
            /* A software graphics driver, as used on a test machine, warns that it is slow. That is not a fault of the page. */
            const driverNote = message.type() === 'warning' && message.text().includes('GL Driver Message');
            if ((message.type() === 'error' || message.type() === 'warning') && !driverNote) {
                problems.push(`console ${message.type()}: ${message.text()}`);
            }
        });
        page.on('pageerror', error => problems.push(`page error: ${error.message}`));
        page.on('requestfailed', request => problems.push(`request failed: ${request.url()}`));
        page.on('response', response => {
            if (response.status() >= 400) {
                problems.push(`http ${response.status()}: ${response.url()}`);
            }
        });
        await page.goto(url, {waitUntil: 'load'});

        /* The numbers written into the text, before the page's scripts replace them. */
        const written = await page.evaluate(async () => {
            const html = await (await fetch('index.html')).text();
            const source = new DOMParser().parseFromString(html, 'text/html');
            return [...source.querySelectorAll('.live')].map(node => [node.dataset.live, node.textContent]);
        });

        /* Visit every figure, so that each one sets itself up and draws. */
        const figures = await page.locator('figure[data-figure]').evaluateAll(nodes => nodes.map(node => node.id));
        if (figures.length < 10) {
            problems.push(`only ${figures.length} figures found`);
        }
        for (const id of figures) {
            await page.locator(`#${id}`).scrollIntoViewIfNeeded();
            await page.waitForTimeout(700);
        }
        /* The slowest figure runs the estimator over a second recording. */
        await page.waitForTimeout(3000);
        for (const id of figures) {
            await page.locator(`#${id}`).scrollIntoViewIfNeeded();
            await page.waitForTimeout(150);
            const drawn = await page.locator(`#${id} canvas`).evaluate(canvas => {
                const context = canvas.getContext('2d');
                const pixels = context.getImageData(0, 0, canvas.width, canvas.height).data;
                let painted = 0;
                for (let i = 3; i < pixels.length; i += 4 * 97) {
                    painted += pixels[i] > 0 ? 1 : 0;
                }
                return {painted, text: canvas.closest('figure').innerText};
            });
            if (drawn.painted < 50) {
                problems.push(`figure ${id} drew almost nothing`);
            }
            if (/Loading|Working|needs WebGL/.test(drawn.text)) {
                problems.push(`figure ${id} did not finish loading`);
            }
        }

        /* Every number in the text must be the estimator's own. */
        const shown = await page.locator('.live').evaluateAll(nodes => nodes.map(node => [node.dataset.live, node.textContent]));
        written.forEach(([name, value], i) => {
            if (shown[i][1] !== value) {
                problems.push(`the text says ${name} = ${value}, but the estimator gives ${shown[i][1]}: update index.html`);
            }
        });

        /* Every phrase which changes a figure must work. */
        const links = await page.locator('button.act').count();
        for (let i = 0; i < links; i++) {
            await page.locator('button.act').nth(i).click();
            await page.waitForTimeout(250);
        }
        await page.waitForTimeout(2500);
        console.log(`${figures.length} figures, ${written.length} live numbers, ${links} links in the text`);
    } finally {
        await browser.close();
        server.close();
    }
    if (problems.length > 0) {
        console.error(problems.join('\n'));
        process.exit(1);
    }
    console.log('The page works.');
}

main().catch(error => {
    console.error(error);
    process.exit(1);
});
